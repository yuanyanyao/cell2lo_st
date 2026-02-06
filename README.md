<<<<<<< HEAD
# st_a-syn_SNpc_C57

Publication-ready organization of spatial transcriptomics analysis for the st_a-syn_SNpc_C57 project.
The primary analysis target is the **bin100 Cell2location** pipeline.

## Repository Structure
- `notebooks/`: Analysis notebooks (`ana.ipynb`, `gem.ipynb`, `model.ipynb`, `test.ipynb`).
- `src/`: Shared Python modules (`config.py`, `data_processing.py`, `plot_utils.py`, `utils.py`).
- `data/`:
  - `data/raw/`: Raw inputs and reference resources.
  - `data/processed/`: Processed inputs (e.g., `bin*_raw.h5ad`, `st.h5ad`).
  - `data/derived/bin100/`: **Primary** derived outputs for bin100 Cell2location analysis.
- `results/`: Model outputs and intermediate analyses.
- `figures/`: All figures (`.png`, `.pdf`).
- `tables/`: All tabular outputs (`.csv`, `.xlsx`).
- `docs/`: Supporting documentation and logs.
- `scripts/`: Utility scripts.
- `metadata/`: Metadata or supplementary files.

## Primary Target: bin100
- Input: `data/processed/bin100_raw.h5ad`
- Derived outputs: `data/derived/bin100/`

## Usage (High-Level)
1. Convert GEF to h5ad using `notebooks/gem.ipynb`.
2. Train Cell2location models with `notebooks/model.ipynb`.
3. Run repair/merge utilities in `notebooks/test.ipynb`.
4. Generate plots and analyses in `notebooks/ana.ipynb`.

## Git LFS
Large `*.h5ad` files are tracked with Git LFS via `.gitattributes`.

## Submission Requirements
I attempted to load the submission requirements from GitHub, but the web tool could not access the repo.
Please paste the requirements here and I will update this README accordingly.
=======
# Stereo-seq CellBin PD analysis

This repo is the working record of my analysis on Stereo-seq **CellBin** spatial transcriptomics data from a unilateral α-syn PD mouse model. The goal is to keep the workflow runnable end-to-end and easy to audit (what was done, with which parameters, and where the figures/tables came from).

A note on terminology: CellBin units are not true single cells. When I write “DA neuron”, “microglia”, etc., I mean **CellBin units enriched for the corresponding marker signals and consistent spatial context**, not a definitive single-cell identity.

---

## What’s here

- `notebooks/Cellbin_analysis_clean.ipynb`  
  The main notebook. It’s cleaned for sharing:
  - English-only comments/section headers
  - no hard-coded personal file paths
  - outputs cleared (so it stays small on GitHub)

- `env/environment.yml` and `requirements.txt`  
  Environment specs used to run the notebook.

- `src/stereo_pd_cellbin/`  
  A small helper module (mostly stats helpers such as Fisher/MWU wrappers and p-value → stars).

- `scripts/run_pipeline.py`  
  A minimal command-line entrypoint to sanity-check loading and run lightweight steps without opening Jupyter.

---

>>>>>>> 13fab32026f08d48b85cd9955c537d0ad280993d
