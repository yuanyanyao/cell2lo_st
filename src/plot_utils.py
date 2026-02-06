"""Plotting utilities to standardize visualization code."""

from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as patches

try:
    import scanpy as sc
except Exception:  # pragma: no cover - optional in some environments
    sc = None


def use_dark_background() -> None:
    """Use dark background for consistent plot styling."""
    plt.style.use("dark_background")


def save_figure(path: str, dpi: int = 300, tight: bool = True) -> None:
    """Save current figure with common options."""
    if tight:
        plt.tight_layout()
    plt.savefig(path, dpi=dpi)


def spatial_plot(
    adata,
    color: str,
    ax=None,
    title: Optional[str] = None,
    cmap: Optional[str] = None,
    spot_size: int = 100,
    vmax: Optional[str] = "p99",
    frameon: bool = False,
) -> None:
    """Wrapper around scanpy spatial plotting with defaults."""
    if sc is None:
        raise ImportError("scanpy is required for spatial_plot")
    sc.pl.spatial(
        adata,
        color=color,
        cmap=cmap,
        spot_size=spot_size,
        vmax=vmax,
        ax=ax,
        show=False,
        title=title,
        frameon=frameon,
    )


def add_roi_boxes(ax, boxes: Sequence[Tuple[int, int, int, int]], edgecolor: str = "cyan") -> None:
    """Add rectangular ROI boxes to an axis.

    Each box is (x_min, x_max, y_min, y_max).
    """
    for x_min, x_max, y_min, y_max in boxes:
        rect = patches.Rectangle(
            (x_min, y_min),
            x_max - x_min,
            y_max - y_min,
            linewidth=2,
            edgecolor=edgecolor,
            facecolor="none",
        )
        ax.add_patch(rect)
