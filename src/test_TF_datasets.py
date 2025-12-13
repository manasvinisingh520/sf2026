"""
Test TF Dataset Overlap Across Brain Regions

This script takes a TF list file and calculates the percent overlap
with genes in each brain region's expression matrix.

Supported TF File Formats:
--------------------------
1. TXT file (recommended for pySCENIC):
   - One TF gene symbol per line

2. CSV file:
   - Tabular format with TFs in a column
   - Use --tf-column to specify column name

3. TSV file:
   - Tab-separated format (same as CSV)

4. Excel file (.xlsx, .xls):
   - Use --tf-column to specify column name
   - Use --no-header if file has no header row

"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd


class TFOverlapAnalyzer:
    """
    Analyze overlap between TF lists and gene expression matrices from different brain regions.
    """
    
    # Available brain regions
    REGIONS = ["EC", "ITG", "PFC", "V1", "V2"]
    
    def __init__(self, data_dir: str = "data", date_prefix: str = "2025-11-16", verbose: bool = True):
        """
        Initialize TF Overlap Analyzer.
        
        Parameters:
        -----------
        data_dir : str, default="data"
            Directory containing matrix files
        date_prefix : str, default="2025-11-16"
            Date prefix for matrix files
        verbose : bool, default=True
            Print progress messages
        """
        self.data_dir = Path(data_dir)
        self.date_prefix = date_prefix
        self.verbose = verbose
    
    def load_tf_list(
        self,
        tf_file: str,
        column: Optional[str] = None,
        header: bool = True
    ) -> Set[str]:
        """
        Load transcription factor list from file.
        
        Supports multiple formats:
        - TXT: One TF per line (default for pySCENIC)
        - CSV: Tabular format with TFs in a column
        - TSV: Tab-separated format
        - Excel: .xlsx or .xls files
        
        Parameters:
        -----------
        tf_file : str
            Path to TF list file
        column : str, optional
            Column name or index to use (for CSV/Excel files).
            If None, uses first column.
        header : bool, default=True
            Whether file has a header row (for CSV/Excel)
        
        Returns:
        --------
        tf_set : set
            Set of TF gene symbols (normalized)
        """
        if not os.path.exists(tf_file):
            raise FileNotFoundError(f"TF file not found: {tf_file}")
        
        if self.verbose:
            print(f"Loading TF list from {tf_file}...")
            print(f"  File format: {Path(tf_file).suffix}")
        
        tfs = set()
        file_ext = Path(tf_file).suffix.lower()
        
        try:
            if file_ext in ['.txt', '.tsv']:
                # Try to detect if it's tab-separated or plain text
                # Check first line for tabs
                with open(tf_file, 'r', encoding='utf-8') as f:
                    first_line = f.readline()
                    has_tabs = '\t' in first_line
                    f.seek(0)  # Reset to beginning of file to reread
                    
                    if has_tabs or column:
                        # Tab-separated file with columns (like TFs_utoronto.txt)
                        df = pd.read_csv(tf_file, sep='\t', header=0 if header else None)
                        
                        if column:
                            if column in df.columns:
                                tf_column = df[column]
                            else:
                                raise ValueError(f"Column '{column}' not found. Available columns: {list(df.columns)}")
                        else:
                            # Try common column names
                            possible_columns = ['HGNC symbol', 'HGNC', 'Gene', 'Gene Symbol', 'Symbol', 'TF']
                            tf_column = None
                            for col_name in possible_columns:
                                if col_name in df.columns:
                                    tf_column = df[col_name]
                                    if self.verbose:
                                        print(f"  Using column: {col_name}")
                                    break
                            
                            if tf_column is None:
                                # Use second column (first is usually index)
                                if len(df.columns) > 1:
                                    tf_column = df.iloc[:, 1]
                                    if self.verbose:
                                        print(f"  Using second column (first column appears to be index)")
                                else:
                                    tf_column = df.iloc[:, 0]
                        
                        for tf in tf_column:
                            if pd.notna(tf):
                                tf_str = str(tf).strip().strip('"').strip("'").strip()
                                if tf_str:
                                    tfs.add(tf_str.upper())
                    else:
                        # Plain text file (one TF per line)
                        # This is the standard pySCENIC format
                        for line in f:
                            tf = line.strip().strip('"').strip("'").strip()
                            if tf and not tf.startswith('#'):  # Skip empty lines and comments
                                tfs.add(tf.upper())  # Normalize to uppercase
            
            elif file_ext == '.csv':
                # CSV file
                df = pd.read_csv(tf_file, header=0 if header else None)
                
                if column:
                    if column in df.columns:
                        tf_column = df[column]
                    else:
                        raise ValueError(f"Column '{column}' not found. Available columns: {list(df.columns)}")
                else:
                    # Use first column
                    tf_column = df.iloc[:, 0]
                
                for tf in tf_column:
                    if pd.notna(tf):
                        tf_str = str(tf).strip().strip('"').strip("'").strip()
                        if tf_str:
                            tfs.add(tf_str.upper())
            
            elif file_ext in ['.xlsx', '.xls']:
                # Excel file
                df = pd.read_excel(tf_file, header=0 if header else None)
                
                if column:
                    if column in df.columns:
                        tf_column = df[column]
                    else:
                        raise ValueError(f"Column '{column}' not found. Available columns: {list(df.columns)}")
                else:
                    # Use first column
                    tf_column = df.iloc[:, 0]
                
                for tf in tf_column:
                    if pd.notna(tf):
                        tf_str = str(tf).strip().strip('"').strip("'").strip()
                        if tf_str:
                            tfs.add(tf_str.upper())
            
            else:
                # Try as plain text as fallback
                if self.verbose:
                    print(f"  Warning: Unknown file extension '{file_ext}', trying as plain text...")
                with open(tf_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        tf = line.strip().strip('"').strip("'").strip()
                        if tf and not tf.startswith('#'):
                            tfs.add(tf.upper())
        
        except Exception as e:
            raise ValueError(f"Error reading TF file '{tf_file}': {e}")
        
        if self.verbose:
            print(f"  Loaded {len(tfs)} unique TFs")
            if len(tfs) > 0:
                sample = sorted(list(tfs))[:5]
                print(f"  Sample TFs: {', '.join(sample)}")
                if len(tfs) > 5:
                    print(f"  ... and {len(tfs) - 5} more")
        
        if len(tfs) == 0:
            raise ValueError(f"No TFs found in file '{tf_file}'. Please check the file format.")
        
        return tfs
    
    def load_region_genes(self, region: str) -> Tuple[Set[str], int, int]:
        """
        Load gene names from a brain region's row annotation file.
        
        Parameters:
        -----------
        region : str
            Brain region (EC, ITG, PFC, V1, V2)
        
        Returns:
        --------
        genes : set
            Set of gene symbols (normalized)
        n_genes : int
            Total number of genes
        n_unique : int
            Number of unique genes
        """
        row_file = self.data_dir / f"{self.date_prefix}_Astrocytes_{region}_row_annotation.txt"
        
        if not row_file.exists():
            # Try alternative date
            alt_file = self.data_dir / f"2025-10-22_Astrocytes_{region}_row_annotation.txt"
            if alt_file.exists():
                row_file = alt_file
            else:
                raise FileNotFoundError(
                    f"Row annotation file not found for region {region}. "
                    f"Tried: {row_file.name} and {alt_file.name}"
                )
        
        if self.verbose:
            print(f"Loading genes from {row_file.name}...")
        
        genes = []
        with open(row_file, 'r') as f:
            for line in f:
                gene = line.strip().strip('"').strip("'").strip()
                if gene:
                    genes.append(gene.upper())  # Normalize to uppercase
        
        gene_set = set(genes)
        n_genes = len(genes)
        n_unique = len(gene_set)
        
        if self.verbose:
            print(f"  Total genes: {n_genes}")
            print(f"  Unique genes: {n_unique}")
            if n_genes != n_unique:
                print(f"  Warning: {n_genes - n_unique} duplicate gene names")
        
        return gene_set, n_genes, n_unique
    
    def calculate_overlap(
        self,
        tf_set: Set[str],
        gene_set: Set[str]
    ) -> Dict[str, float]:
        """
        Calculate overlap statistics between TF list and gene set.
        
        Parameters:
        -----------
        tf_set : set
            Set of TF gene symbols
        gene_set : set
            Set of genes from expression matrix
        
        Returns:
        --------
        stats : dict
            Dictionary with overlap statistics
        """
        # Find overlapping genes
        overlap = tf_set.intersection(gene_set)
        
        # Calculate statistics
        n_tfs = len(tf_set)
        n_genes = len(gene_set)
        n_overlap = len(overlap)
        
        if n_tfs == 0:
            percent_tf_in_genes = 0.0
        else:
            percent_tf_in_genes = (n_overlap / n_tfs) * 100
        
        if n_genes == 0:
            percent_genes_are_tf = 0.0
        else:
            percent_genes_are_tf = (n_overlap / n_genes) * 100
        
        stats = {
            "n_tfs": n_tfs,
            "n_genes": n_genes,
            "n_overlap": n_overlap,
            "percent_tf_in_genes": percent_tf_in_genes,
            "percent_genes_are_tf": percent_genes_are_tf,
            "overlap_genes": overlap,
            "missing_tfs": tf_set - gene_set,
        }
        
        return stats
    
    def analyze_all_regions(
        self,
        tf_file: str,
        tf_column: Optional[str] = None,
        tf_header: bool = True
    ) -> pd.DataFrame:
        """
        Analyze TF overlap across all brain regions.
        
        Parameters:
        -----------
        tf_file : str
            Path to TF list file (supports .txt, .csv, .tsv, .xlsx formats)
        tf_column : str, optional
            Column name to use for CSV/Excel files (default: first column)
        tf_header : bool, default=True
            Whether file has a header row (for CSV/Excel)
        
        Returns:
        --------
        results_df : pd.DataFrame
            DataFrame with overlap statistics for each region
        """
        # Load TF list
        tf_set = self.load_tf_list(tf_file, column=tf_column, header=tf_header)
        
        if self.verbose:
            print("\n" + "=" * 60)
            print("Analyzing overlap across brain regions")
            print("=" * 60)
        
        results = []
        detailed_results = {}
        
        # Analyze each region
        for region in self.REGIONS:
            if self.verbose:
                print(f"\nRegion: {region}")
                print("-" * 40)
            
            try:
                # Load genes for this region
                gene_set, n_genes, n_unique = self.load_region_genes(region)
                
                # Calculate overlap
                stats = self.calculate_overlap(tf_set, gene_set)
                
                # Store results
                result_row = {
                    "Region": region,
                    "Total_Genes": n_unique,
                    "Total_TFs": stats["n_tfs"],
                    "Overlapping_TFs": stats["n_overlap"],
                    "Percent_TF_in_Genes": round(stats["percent_tf_in_genes"], 2),
                    "Percent_Genes_are_TF": round(stats["percent_genes_are_tf"], 2),
                    "Missing_TFs": len(stats["missing_tfs"]),
                }
                results.append(result_row)
                detailed_results[region] = stats
                
                if self.verbose:
                    print(f"  Overlapping TFs: {stats['n_overlap']}/{stats['n_tfs']} "
                          f"({stats['percent_tf_in_genes']:.2f}%)")
                    print(f"  TFs in dataset: {stats['n_overlap']}/{n_unique} "
                          f"({stats['percent_genes_are_tf']:.2f}%)")
                    if stats['missing_tfs']:
                        print(f"  Missing TFs: {stats['missing_tfs']}")
                        if len(stats['missing_tfs']) <= 10:
                            print(f"    {', '.join(list(stats['missing_tfs'])[:10])}")
                        else:
                            print(f"    {', '.join(list(stats['missing_tfs'])[:10])} ... "
                                  f"({len(stats['missing_tfs'])} total)")
            
            except FileNotFoundError as e:
                if self.verbose:
                    print(f"  Error: {e}")
                results.append({
                    "Region": region,
                    "Total_Genes": None,
                    "Total_TFs": stats["n_tfs"] if 'stats' in locals() else len(tf_set),
                    "Overlapping_TFs": None,
                    "Percent_TF_in_Genes": None,
                    "Percent_Genes_are_TF": None,
                    "Missing_TFs": None,
                })
        
        # Create DataFrame
        results_df = pd.DataFrame(results)
        
        # Store detailed results for access
        self.detailed_results = detailed_results
        self.tf_set = tf_set
        
        return results_df
    
    def print_summary(self, results_df: pd.DataFrame):
        """
        Print a simple summary with just percentages for each region.
        
        Parameters:
        -----------
        results_df : pd.DataFrame
            Results DataFrame from analyze_all_regions()
        """
        print("\nTF Overlap Percentages by Region:")
        print("-" * 40)
        for _, row in results_df.iterrows():
            if pd.notna(row['Percent_TF_in_Genes']):
                print(f"{row['Region']}: {row['Percent_TF_in_Genes']:.2f}%")
            else:
                print(f"{row['Region']}: N/A")
    
    def save_results(self, results_df: pd.DataFrame, output_file: str = "tf_overlap_results.csv"):
        """
        Save results to CSV file.
        
        Parameters:
        -----------
        results_df : pd.DataFrame
            Results DataFrame
        output_file : str, default="tf_overlap_results.csv"
            Output file path
        """
        results_df.to_csv(output_file, index=False)
        if self.verbose:
            print(f"\nResults saved to: {output_file}")
    
    def get_missing_tfs_by_region(self) -> Dict[str, Set[str]]:
        """
        Get missing TFs for each region.
        
        Returns:
        --------
        missing_by_region : dict
            Dictionary mapping region names to sets of missing TF symbols
        """
        if not hasattr(self, 'detailed_results'):
            raise ValueError("Must run analyze_all_regions() first")
        
        return {
            region: stats["missing_tfs"]
            for region, stats in self.detailed_results.items()
            if "missing_tfs" in stats
        }
    
    def get_overlap_by_region(self) -> Dict[str, Set[str]]:
        """
        Get overlapping TFs for each region.
        
        Returns:
        --------
        overlap_by_region : dict
            Dictionary mapping region names to sets of overlapping TF symbols
        """
        if not hasattr(self, 'detailed_results'):
            raise ValueError("Must run analyze_all_regions() first")
        
        return {
            region: stats["overlap_genes"]
            for region, stats in self.detailed_results.items()
            if "overlap_genes" in stats
        }


def main():
    """Main function to run TF overlap analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description=(
            "Analyze TF overlap with gene expression matrices from different brain regions.\n"
            "Supports .txt (one TF per line - standard pySCENIC format), .csv, .tsv, and .xlsx files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--tf-file",
        type=str,
        help="Path to TF list file. Supports: .txt (one per line), .csv, .tsv, .xlsx"
    )
    parser.add_argument(
        "--tf-column",
        type=str,
        default=None,
        help="Column name/index to use for CSV/Excel files (default: first column)"
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="File has no header row (for CSV/Excel)"
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="data",
        help="Directory containing matrix files (default: data)"
    )
    parser.add_argument(
        "--date",
        type=str,
        default="2025-11-16",
        help="Date prefix for matrix files (default: 2025-11-16)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional: Output CSV file to save results"
    )
    parser.add_argument(
        "--regions",
        type=str,
        nargs="+",
        choices=["EC", "ITG", "PFC", "V1", "V2"],
        help="Specific regions to analyze (default: all)"
    )
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = TFOverlapAnalyzer(
        data_dir=args.data_dir,
        date_prefix=args.date,
        verbose=False  # Minimal output by default
    )
    
    # Override regions if specified
    if args.regions:
        analyzer.REGIONS = args.regions
    
    # Run analysis
    results_df = analyzer.analyze_all_regions(
        args.tf_file,
        tf_column=args.tf_column,
        tf_header=not args.no_header
    )
    
    # Print summary (just percentages)
    analyzer.print_summary(results_df)
    
    # Save results (optional)
    if args.output:
        analyzer.save_results(results_df, args.output)


if __name__ == "__main__":
    main()

