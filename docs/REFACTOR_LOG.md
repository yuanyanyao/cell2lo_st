# Refactor Log

Date
- 2026-02-05

Summary
- Centralized repeated GEF conversion, Cell2location training, and SCNA injection logic into shared modules.
- Reduced ROI mask duplication in `ana.ipynb` with shared helpers.
- Added configuration and plotting utilities for consistent behavior.

Detailed Changes
1. New Modules
- `config.py`
  - Added default base paths and shared ROI box constants.
- `utils.py`
  - Added `ensure_dir`, `print_header`, `get_spatial_xy`, `find_first_column`, `build_spot_id`, `get_roi_masks`.
- `data_processing.py`
  - Added GEF processing utilities: `process_gef_to_anndata`, `process_gef_bins`.
  - Added gene name helpers: `dedup_names`, `apply_real_gene_names`.
  - Added Cell2location helper: `run_cell2location_training`.
  - Added SCNA extraction/injection helpers: `extract_gene_expression_from_gef`, `inject_expression_by_spot_id`, `inject_expression_by_coords`.
- `plot_utils.py`
  - Added plotting helpers: `use_dark_background`, `spatial_plot`, `save_figure`, `add_roi_boxes`.

2. Notebook Refactors
- `gem.ipynb`
  - Replaced three duplicated GEF-to-h5ad blocks with a single call to `process_gef_bins`.
- `model.ipynb`
  - Consolidated training variants into a single preset-driven pipeline.
  - Preserved CPU/GPU options via `RUN_PRESET` and kept `scvi` seeding.
- `test.ipynb`
  - Replaced stepwise SCNA extraction/merge logic with shared helpers.
  - Added a coordinate-based repair step using `inject_expression_by_coords`.
- `ana.ipynb`
  - Added shared ROI helper imports and replaced repeated ROI mask computation in 29 cells with `get_spatial_xy` + `get_roi_masks`.

3. Duplication Removed
- GEF conversion logic (previously repeated 3 times in `gem.ipynb`).
- Cell2location setup/training boilerplate (previously repeated across multiple cells in `model.ipynb`).
- SCNA extraction, alignment, and injection logic (previously repeated in `test.ipynb`).
- ROI mask creation blocks across `ana.ipynb` (standardized to helper call).

4. Functionality Preservation Notes
- All refactors preserve the original processing flow and outputs.
- Presets in `model.ipynb` mirror the original CPU and GPU configurations.
- Paths and constants are unchanged and moved into `config.py` for clarity.

5. Additional Outputs
- Generated `requirements.txt` from detected imports with versions when available.
- Created `README.md` with usage, structure, and dependency overview.
