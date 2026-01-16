"""
Script to extract 'HGNC symbol' column from TFs_utoronto.csv and save as text file.
Only saves TFs that appear in the region's matrix file.
Outputs one gene symbol per line.
"""

import pandas as pd
import os
import argparse

def convert_tfs_csv_to_txt(region, cell_type="Astrocytes"):
    """Extract 'HGNC symbol' column from TFs_utoronto.csv and save as text file.
    Only includes TFs that appear in the region's matrix file.
    
    Parameters:
    -----------
    region : str
        Region name (EC, ITG, PFC, V1, V2, CrossRegion)
    cell_type : str
        Cell type (e.g., 'Astrocytes', 'Microglia'). Default: 'Astrocytes'
    """
    
    input_file = "data/TFs_utoronto.csv"
    
    # Set output directory based on cell_type
    output_dir = os.path.join("data", "GRN", cell_type)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get matrix filename based on cell_type (same pattern as export_matrices_to_csv.py)
    if cell_type == "Microglia":
        matrix_file = os.path.join(output_dir, f"{region}_Microglia_matrix.csv")
    else:
        matrix_file = os.path.join(output_dir, f"{region}_matrix.csv")
    
    output_file = os.path.join(output_dir, f"{region}_TFs.txt")
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        return
    
    # Check if matrix file exists
    if not os.path.exists(matrix_file):
        print(f"Error: Matrix file not found: {matrix_file}")
        return
    
    print(f"Reading matrix file: {matrix_file}")
    # Read V2_matrix.csv and extract gene names (first column, no header)
    df = pd.read_csv(matrix_file)
    matrix_genes = set(df.iloc[:, 0])
    
    print(f"Found {len(matrix_genes)} genes in {matrix_file}")
    
    print(f"\nReading CSV file: {input_file}")
    
    # Read CSV file
    df = pd.read_csv(input_file)
    
    print(f"Loaded {len(df)} rows and {len(df.columns)} columns")
    
    # Check if 'HGNC symbol' column exists
    if 'HGNC symbol' not in df.columns:
        print(f"Error: 'HGNC symbol' column not found. Available columns: {list(df.columns)}")
        return
    
    # Extract only the 'HGNC symbol' column and remove NaN values
    hgnc_symbols = df['HGNC symbol'].dropna().unique()
    
    print(f"Found {len(hgnc_symbols)} unique HGNC symbols in TFs file")
    
    # Filter to only include TFs that appear in V2_matrix.csv
    filtered_symbols = [symbol for symbol in hgnc_symbols if symbol in matrix_genes]
    
    print(f"\nFiltered to {len(filtered_symbols)} TFs that appear in {os.path.basename(matrix_file)}")
    print(f"  (Removed {len(hgnc_symbols) - len(filtered_symbols)} TFs not in matrix)")
    
    # Save as text file (one symbol per line)
    print(f"\nWriting to: {output_file}")
    with open(output_file, 'w') as f:
        for symbol in sorted(filtered_symbols):  # Sort for consistency
            f.write(f"{symbol}\n")
    
    # Check file size
    file_size = os.path.getsize(output_file) / 1024  # KB
    print(f"  Successfully saved filtered HGNC symbols")
    print(f"  Output file: {output_file}")
    print(f"  File size: {file_size:.2f} KB")
    print(f"  Number of symbols: {len(filtered_symbols)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Extract TFs from TFs_utoronto.csv that appear in region matrix file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '-region',
        type=str,
        default='EC',
        help='Region to process (e.g., EC, ITG, PFC, V1, V2, CrossRegion). Default: EC'
    )
    parser.add_argument(
        '-cell_type',
        type=str,
        default="Astrocytes",
        help='Cell type (e.g., "Astrocytes", "Microglia"). Default: "Astrocytes"'
    )
    
    args = parser.parse_args()
    
    print(f"Processing region: {args.region}")
    print(f"Cell type: {args.cell_type}")
    print("="*60)
    
    convert_tfs_csv_to_txt(region=args.region, cell_type=args.cell_type)

