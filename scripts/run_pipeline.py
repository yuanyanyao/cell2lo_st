"""
Run the cleaned CellBin analysis in a non-interactive way.

This script is intentionally minimal: it loads an AnnData object and
writes a copy to results/ as a sanity check. For full analysis, use the notebook.

Usage:
  python scripts/run_pipeline.py --input data/raw/your_file.h5ad --outdir results
"""
from __future__ import annotations

import argparse
import os
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to input .h5ad")
    ap.add_argument("--outdir", default="results", help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    adata = sc.read_h5ad(args.input)
    print(adata)

    out_path = os.path.join(args.outdir, "adata_loaded.h5ad")
    adata.write(out_path)
    print(f"[OK] Wrote: {out_path}")


if __name__ == "__main__":
    main()
