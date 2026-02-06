"""General utility helpers used across notebooks and scripts."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np


def ensure_dir(path: Path | str) -> Path:
    """Create a directory if it does not exist and return it as Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def print_header(title: str, char: str = "=") -> None:
    """Print a formatted header for readability in notebooks."""
    line = char * 80
    print(f"\n{line}\n{title}\n{line}")


def get_spatial_xy(adata) -> Tuple[np.ndarray, np.ndarray]:
    """Return x, y spatial coordinates from AnnData."""
    coords = adata.obsm.get("spatial")
    if coords is None:
        raise ValueError("AnnData is missing obsm['spatial']")
    return coords[:, 0], coords[:, 1]


def find_first_column(columns: Sequence[str], substrings: Iterable[str]) -> Optional[str]:
    """Find the first column name containing any of the substrings."""
    for col in columns:
        for sub in substrings:
            if sub in col:
                return col
    return None


def build_spot_id(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Build spot_id strings from x and y arrays."""
    return x.astype(int).astype(str) + "_" + y.astype(int).astype(str)


def get_roi_masks(
    x: np.ndarray,
    y: np.ndarray,
    left_box: Tuple[int, int, int, int],
    right_box: Tuple[int, int, int, int],
) -> Tuple[np.ndarray, np.ndarray]:
    """Return boolean masks for left/right ROI boxes.

    Boxes are (x_min, x_max, y_min, y_max).
    """
    x_min_l, x_max_l, y_min_l, y_max_l = left_box
    x_min_r, x_max_r, y_min_r, y_max_r = right_box

    mask_l = (x >= x_min_l) & (x <= x_max_l) & (y >= y_min_l) & (y <= y_max_l)
    mask_r = (x >= x_min_r) & (x <= x_max_r) & (y >= y_min_r) & (y <= y_max_r)
    return mask_l, mask_r
