"""Data processing utilities shared across notebooks."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse

from utils import ensure_dir

try:
    import stereo as st
except Exception:  # pragma: no cover
    st = None


def rotate_coords_90_clockwise(coords: np.ndarray) -> np.ndarray:
    """Rotate (x, y) -> (y, -x) and shift to non-negative."""
    new_coords = np.zeros_like(coords)
    new_coords[:, 0] = coords[:, 1]
    new_coords[:, 1] = -coords[:, 0]
    new_coords[:, 0] -= new_coords[:, 0].min()
    new_coords[:, 1] -= new_coords[:, 1].min()
    return new_coords


def dedup_names(names: Iterable[str]) -> List[str]:
    """Deduplicate gene names by appending -N suffixes."""
    seen: Dict[str, int] = {}
    out: List[str] = []
    for name in map(str, names):
        if name not in seen:
            seen[name] = 0
            out.append(name)
        else:
            seen[name] += 1
            out.append(f"{name}-{seen[name]}")
    return out


def apply_real_gene_names(adata, real_gene_names: Optional[Iterable[str]]) -> None:
    """Replace adata.var_names with real gene names if lengths match."""
    if real_gene_names is None:
        return
    real_gene_names = list(real_gene_names)
    if len(real_gene_names) == adata.shape[1]:
        adata.var_names = real_gene_names


def process_gef_to_anndata(
    gef_path: str | Path,
    bin_size: int,
    output_path: str | Path,
    min_gene: int = 10,
    rotate: bool = True,
) -> None:
    """Load GEF, rotate coords, filter, convert to AnnData, fix gene names, and save."""
    if st is None:
        raise ImportError("stereo is required for process_gef_to_anndata")

    gef_path = str(gef_path)
    output_path = Path(output_path)
    ensure_dir(output_path.parent)

    print(f"正在加载 GEF 文件 (Bin{bin_size})...")
    data = st.io.read_gef(file_path=gef_path, bin_size=bin_size)
    print(f"原始数据加载完成。Shape: {data.shape}")

    real_gene_names = None
    if hasattr(data.genes, "real_gene_name") and data.genes.real_gene_name is not None:
        print("\n>>> 检测到 real_gene_name，已提取到内存。")
        real_gene_names = data.genes.real_gene_name
    else:
        print("⚠️ 未找到 real_gene_name，稍后将继续使用 ID。")

    if rotate:
        print("\n正在旋转坐标...")
        data.position = rotate_coords_90_clockwise(data.position)

    data.tl.cal_qc()
    data.tl.filter_cells(min_gene=min_gene, inplace=True)
    print(f"过滤后 Cell 数: {data.shape[0]}")

    print("\n正在转换为 AnnData...")
    adata = st.io.stereo_to_anndata(data, flavor="scanpy")

    apply_real_gene_names(adata, real_gene_names)
    print("正在检查并处理重复基因名...")
    adata.var_names = dedup_names(adata.var_names)
    print(f"✅ 最终基因名预览: {adata.var_names[:5].tolist()}")

    if not scipy.sparse.issparse(adata.X):
        adata.X = np.round(adata.X)

    adata.layers["counts"] = adata.X.copy()

    print(f"\n正在保存至: {output_path}")
    adata.write(output_path)
    print("🎉 成功！")


def process_gef_bins(
    gef_path: str | Path,
    bin_sizes: Sequence[int],
    output_dir: str | Path,
    min_gene: int = 10,
    rotate: bool = True,
) -> None:
    """Process multiple GEF bin sizes to h5ad outputs."""
    output_dir = Path(output_dir)
    ensure_dir(output_dir)
    for bin_size in bin_sizes:
        output_path = output_dir / f"bin{bin_size}_raw.h5ad"
        process_gef_to_anndata(
            gef_path=gef_path,
            bin_size=bin_size,
            output_path=output_path,
            min_gene=min_gene,
            rotate=rotate,
        )


def prepare_c2l_inputs(adata, ref_df: pd.DataFrame) -> Tuple[object, pd.DataFrame]:
    """Intersect genes between adata and reference and return aligned copies."""
    intersect_genes = np.intersect1d(adata.var_names, ref_df.index)
    adata_c2l = adata[:, intersect_genes].copy()
    current_ref = ref_df.loc[intersect_genes, :].copy()
    return adata_c2l, current_ref


def run_cell2location_training(
    adata,
    ref_df: pd.DataFrame,
    n_cells_mean: int,
    max_epochs: int,
    batch_size: Optional[int] = None,
    train_size: float = 1,
    detection_alpha: int = 20,
    accelerator: Optional[str] = None,
    save_dir: Optional[str | Path] = None,
    export_posterior: bool = False,
    export_num_samples: int = 1000,
    export_batch_size: Optional[int] = None,
    use_gpu_export: Optional[bool] = None,
) -> Tuple[object, Optional[dict]]:
    """Train a Cell2location model and optionally export posterior results.

    Returns (model, posterior_dict_or_None).
    """
    import cell2location

    adata_c2l, current_ref = prepare_c2l_inputs(adata, ref_df)
    cell2location.models.Cell2location.setup_anndata(adata_c2l, batch_key=None)

    mod = cell2location.models.Cell2location(
        adata_c2l,
        cell_state_df=current_ref,
        detection_alpha=detection_alpha,
        N_cells_per_location=n_cells_mean,
    )

    train_kwargs = {
        "max_epochs": max_epochs,
        "batch_size": batch_size,
        "train_size": train_size,
    }
    if accelerator is not None:
        train_kwargs["accelerator"] = accelerator

    mod.train(**train_kwargs)

    posterior = None
    if export_posterior:
        posterior = mod.export_posterior(
            num_samples=export_num_samples,
            batch_size=export_batch_size or batch_size,
            use_gpu=use_gpu_export if use_gpu_export is not None else False,
        )

    if save_dir is not None:
        mod.save(str(save_dir), overwrite=True)

    return mod, posterior


def extract_gene_expression_from_gef(
    gef_path: str | Path,
    bin_size: int,
    target_gene: str,
) -> pd.DataFrame:
    """Extract expression values for a target gene from a GEF file."""
    if st is None:
        raise ImportError("stereo is required for extract_gene_expression_from_gef")

    data = st.io.read_gef(file_path=str(gef_path), bin_size=bin_size)
    if data.genes.real_gene_name is not None:
        all_genes = data.genes.real_gene_name
    else:
        all_genes = data.gene_names

    found_gene = None
    for g in all_genes:
        if target_gene in str(g):
            found_gene = str(g)
            break

    if found_gene is None:
        raise ValueError(f"Gene '{target_gene}' not found in GEF")

    if data.genes.real_gene_name is not None:
        gene_idx = np.where(data.genes.real_gene_name == found_gene)[0][0]
    else:
        gene_idx = np.where(data.gene_names == found_gene)[0][0]

    if scipy.sparse.issparse(data.exp_matrix):
        expr_values = data.exp_matrix[:, gene_idx].toarray().flatten()
    else:
        expr_values = data.exp_matrix[:, gene_idx].flatten()

    df = pd.DataFrame(
        {
            "x": data.position[:, 0],
            "y": data.position[:, 1],
            "gene_name": found_gene,
            "expression": expr_values,
        }
    )
    df["spot_id"] = df["x"].astype(int).astype(str) + "_" + df["y"].astype(int).astype(str)
    return df


def inject_expression_by_spot_id(
    adata,
    expr_df: pd.DataFrame,
    target_obs_key: str,
) -> np.ndarray:
    """Inject expression into adata.obs by spot_id alignment."""
    expr_df = expr_df.copy()
    expr_df["spot_id"] = expr_df["x"].astype(int).astype(str) + "_" + expr_df["y"].astype(int).astype(str)
    expr_map = dict(zip(expr_df["spot_id"], expr_df["expression"]))

    new_expr = np.zeros(adata.n_obs)
    match_count = 0
    for i, idx in enumerate(adata.obs_names):
        clean_idx = str(idx).split("-")[0].replace(":", "_")
        if clean_idx in expr_map:
            new_expr[i] = expr_map[clean_idx]
            match_count += 1
    adata.obs[target_obs_key] = new_expr
    return new_expr


def inject_expression_by_coords(
    adata,
    expr_df: pd.DataFrame,
    target_obs_key: str,
    use_rotation: bool = False,
    threshold_small: float = 2.0,
    threshold_large: float = 150.0,
) -> np.ndarray:
    """Inject expression into adata.obs by coordinate matching with KDTree."""
    from scipy.spatial import cKDTree

    coords_target = adata.obsm.get("spatial")
    if coords_target is None:
        raise ValueError("AnnData missing obsm['spatial']")

    coords_source = expr_df[["x", "y"]].values.astype(float)

    if use_rotation:
        coords_source = rotate_coords_90_clockwise(coords_source)

    tree = cKDTree(coords_source)
    distances, indices = tree.query(coords_target, k=1)

    coord_span = coords_target[:, 0].max() - coords_target[:, 0].min()
    threshold = threshold_large if coord_span > 2000 else threshold_small

    matched_mask = distances < threshold

    new_expr = np.zeros(adata.n_obs)
    valid_indices = indices[matched_mask]
    valid_values = expr_df["expression"].values[valid_indices]
    new_expr[matched_mask] = valid_values

    adata.obs[target_obs_key] = new_expr
    return new_expr
