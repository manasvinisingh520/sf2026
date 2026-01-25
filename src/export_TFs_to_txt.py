"""Extract 'HGNC symbol' from TFs_utoronto.csv and save as text file. Only includes TFs in region matrix."""

import pandas as pd
import os
import argparse

def convert_tfs_csv_to_txt(region, cell_type="Astrocytes"):
    """Extract TFs that appear in the region's matrix file and save to text."""
    input_file = "data/TFs_utoronto.csv"
    output_dir = os.path.join("data", "GRN", cell_type)
    os.makedirs(output_dir, exist_ok=True)
    
    # Get file paths
    is_microglia = cell_type == "Microglia"
    matrix_file = os.path.join(output_dir, f"{region}_Microglia_matrix.csv" if is_microglia else f"{region}_matrix.csv")
    output_file = os.path.join(output_dir, f"{region}_TFs.txt")
    
    # Check files exist
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        return
    if not os.path.exists(matrix_file):
        print(f"Error: Matrix file not found: {matrix_file}")
        return
    
    # Read matrix genes (first column, no header)
    matrix_genes = set(pd.read_csv(matrix_file).iloc[:, 0])
    print(f"Found {len(matrix_genes)} genes in matrix")
    
    # Read and filter TFs
    df = pd.read_csv(input_file)
    if 'HGNC symbol' not in df.columns:
        print(f"Error: 'HGNC symbol' column not found. Available: {list(df.columns)}")
        return
    
    hgnc_symbols = df['HGNC symbol'].dropna().unique()
    filtered_symbols = sorted([s for s in hgnc_symbols if s in matrix_genes])
    
    # Save filtered TFs
    with open(output_file, 'w') as f:
        f.write('\n'.join(filtered_symbols) + '\n')
    
    print(f"Saved {len(filtered_symbols)} TFs (removed {len(hgnc_symbols) - len(filtered_symbols)} not in matrix)")
    print(f"Output: {output_file} ({os.path.getsize(output_file) / 1024:.2f} KB)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract TFs from TFs_utoronto.csv that appear in region matrix file')
    parser.add_argument('-region', type=str, default='EC', help='Region to process (default: EC)')
    parser.add_argument('-cell_type', type=str, default="Astrocytes", help='Cell type (default: Astrocytes)')
    
    args = parser.parse_args()
    print(f"Processing {args.region} ({args.cell_type})\n{'='*60}")
    convert_tfs_csv_to_txt(region=args.region, cell_type=args.cell_type)

