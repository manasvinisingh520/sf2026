"""
Configuration constants for single-cell data processing.
"""

# Mapping of region names to Excel sheet tab indices (0-based)
REGION_TO_TAB = {
    "EC": 0,
    "ITG": 1,
    "PFC": 2,
    "V2": 3,
    "V1": 4,
}

# Valid region choices
REGIONS = list(REGION_TO_TAB.keys())

# Default file naming patterns
ASTROCYTES_BASE_PREFIX = "2025-10-22_Astrocytes_{region}"
DEFAULT_METADATA_PATH = "data/2025-11-16_Astrocytes_metadata.xlsx"

# Date prefixes for file naming
METADATA_DATE_PREFIX = "2025-11-16"
MATRIX_DATE_PREFIX = "2025-10-22"

# Directory paths (relative to project root)
DATA_DIR = "data"
DGE_RESULTS_BASE_DIR = "results"

