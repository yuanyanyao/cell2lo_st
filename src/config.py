"""Project configuration and shared constants."""

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
ANALYSIS_RESULTS_DIR = BASE_DIR / "analysis_results"
PROCESSED_DATA_DIR = BASE_DIR / "processed_data"
REFER_DIR = BASE_DIR / "refer"
FIGURES_DIR = BASE_DIR / "figures"

# Default input paths (override in notebooks as needed)
DEFAULT_GEF_PATH = Path(
    "/home/yuanyanyao/Downloads/hsyn/Result_X101SC25093716-Z01-F002/1.SAW/BrPD1/feature_expression/Y01504E1.tissue.gef"
)
DEFAULT_REFERENCE_CSV = REFER_DIR / "Whole_Brain_Signature_Final_Symbols.csv"

# Standard ROI boxes used in many spatial plots (x_min, x_max, y_min, y_max)
ROI_LEFT_BOX = (6000, 7800, 11400, 13200)   # Model (Left)
ROI_RIGHT_BOX = (12200, 14000, 11400, 13200)  # Control (Right)
