"""
Compute s_impact and Z statistics for all DGE files in a given region.

For each DGE file (condition comparison) in the specified region:
- Calculate s_impact = sum of (padj < 0.05) * |log2FoldChange|
- Calculate Z = log2FoldChange / lfcSE for all genes
"""

import pandas as pd
import numpy as np
import argparse
import re
from pathlib import Path
from utils import get_dge_results_dir
from config import REGIONS, DGE_RESULTS_BASE_DIR


def compute_s_impact_and_z(dge_file: Path):
    """
    Compute s_impact and Z statistics for a single DGE file.
    
    Parameters:
    -----------
    dge_file : Path
        Path to DGE results CSV file
        
    Returns:
    --------
    s_impact : float
        Sum of |log2FoldChange| for genes with padj < 0.05
    z_values : np.ndarray
        Array of Z values (log2FoldChange / lfcSE) for all genes
    """
    # Read the DGE results file
    df = pd.read_csv(dge_file, index_col=0)
    
    # Calculate s_impact: sum of |log2FoldChange| for genes with padj < 0.05
    # Filter out NaN padj values
    significant_mask = (df['padj'] < 0.05) & (df['padj'].notna())
    s_impact = (df.loc[significant_mask, 'log2FoldChange'].abs()).sum()
    
    # Calculate Z = log2FoldChange / lfcSE for all genes
    # Handle division by zero or NaN lfcSE values
    z_values = df['log2FoldChange'] / df['lfcSE'].replace(0, np.nan)
    
    return s_impact, z_values


def main():
    parser = argparse.ArgumentParser(
        description='Compute s_impact and Z statistics for all DGE files in a region',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-version', type=str, choices=['1', '2', '3', '4', 'final'], 
                       default='final', 
                       help='DGE version (1, 2, 3, 4, or final). Default: final')
    parser.add_argument('-region', type=str, choices=REGIONS + ['CrossRegion'], required=True, 
                       help='Region name (e.g., EC, ITG, PFC, V1, V2)')
    
    args = parser.parse_args()
    
    # Get results directory
    if args.version == 'final':
        results_dir = Path(DGE_RESULTS_BASE_DIR) / 'dge_final'
    else:
        version = int(args.version)
        results_dir = get_dge_results_dir(version)
    
    # Find all DGE files for the specified region
    # Pattern matches files like: dge_results_EC_*.csv or dge_results_microglia_EC_*.csv
    pattern1 = f"dge_results_microglia_{args.region}_*.csv"
    
    matching_files = sorted(list(results_dir.glob(pattern1)))
    
    if len(matching_files) == 0:
        raise FileNotFoundError(
            f"No DGE files found for region '{args.region}'\n"
            f"Searched in directory: {results_dir.absolute()}\n"
            f"Patterns: {pattern1}"
        )
    
    print(f"Found {len(matching_files)} DGE files for region '{args.region}'")
    print(f"Results directory: {results_dir.absolute()}\n")
    
    # Initialize lists to store results
    s_impact_list = []
    z_list = []  # List of arrays, one per file
    distances = []  # List of distances between conditions
    
    # Process each file
    for dge_file in matching_files:
        # Extract condition numbers from filename (e.g., "1_vs_2" from "dge_results_EC_..._1_vs_2.csv")
        match = re.search(r'(\d+)_vs_(\d+)', dge_file.name)
        if match:
            cond1 = int(match.group(1))
            cond2 = int(match.group(2))
            distance = abs(cond1 - cond2)
        else:
            # If pattern not found, default to distance 1
            distance = 1
            print(f"Warning: Could not extract condition numbers from {dge_file.name}, using distance=1")
        
        distances.append(distance)
        s_impact, z_values = compute_s_impact_and_z(dge_file)
        s_impact_list.append(s_impact)
        z_list.append(z_values.values)  # Convert to numpy array
        
    # Compute weighted averages
    distances = np.array(distances)
    s_impact_array = np.array(s_impact_list)
    z_medians_array = np.array([np.nanmedian(np.abs(z)) for z in z_list])
    
    # total_s_impact = sum(d(c) * s_impact(c)) / sum(d(c))
    total_s_impact = np.sum(distances * s_impact_array) / np.sum(distances)
    
    # total_z = sum(d(c) * z_median(c)) / sum(d(c))
    total_z = np.sum(distances * z_medians_array) / np.sum(distances)
        
    # Print both lists
    print("=" * 60)
    print("Summary Results")
    print("=" * 60)
    print("\ns_impact list:")
    
    print("\nFull s_impact list (for copy-paste):")
    # Convert numpy types to native Python types
    s_impact_list_python = [float(x) for x in s_impact_list]
    print(s_impact_list_python)
    
    print("\nZ list summary (median(|Z|) per file):")
    z_medians = [float(np.nanmedian(np.abs(z))) for z in z_list]
    print(z_medians)
    
    print("\n" + "=" * 60)
    print("Weighted Averages (weighted by condition distance)")
    print("=" * 60)
    print(f"\ntotal_s_impact = {float(total_s_impact):.4f}")
    print(f"total_z (median(|Z|)) = {float(total_z):.4f}")


if __name__ == "__main__":
    main()
