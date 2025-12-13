"""
Test Motif Dataset Overlap with TF List

This script tests the overlap between a TF list and motif databases,
checking which TFs have corresponding motifs in the database.
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd


class MotifOverlapAnalyzer:
    """
    Analyze overlap between TF lists and motif databases.
    """
    
    def __init__(self, verbose: bool = True):
        """
        Initialize Motif Overlap Analyzer.
        
        Parameters:
        -----------
        verbose : bool, default=True
            Print progress messages
        """
        self.verbose = verbose
    
    def load_tf_list(
        self,
        tf_file: str,
        column: Optional[str] = None,
        header: bool = True
    ) -> Set[str]:
        """
        Load transcription factor list from file.
        
        Supports multiple formats: .txt, .csv, .tsv, .xlsx
        
        Parameters:
        -----------
        tf_file : str
            Path to TF list file
        column : str, optional
            Column name to use for CSV/Excel files
        header : bool, default=True
            Whether file has a header row
        
        Returns:
        --------
        tf_set : set
            Set of TF gene symbols (normalized to uppercase)
        """
        if not os.path.exists(tf_file):
            raise FileNotFoundError(f"TF file not found: {tf_file}")
        
        if self.verbose:
            print(f"Loading TF list from {tf_file}...")
        
        tfs = set()
        file_ext = Path(tf_file).suffix.lower()
        
        try:
            if file_ext in ['.txt', '.tsv']:
                # Check if tab-separated or plain text
                with open(tf_file, 'r', encoding='utf-8') as f:
                    first_line = f.readline()
                    has_tabs = '\t' in first_line
                    f.seek(0)
                    
                    if has_tabs or column:
                        # Tab-separated with columns
                        df = pd.read_csv(tf_file, sep='\t', header=0 if header else None)
                        
                        if column:
                            if column in df.columns:
                                tf_column = df[column]
                            else:
                                raise ValueError(f"Column '{column}' not found. Available: {list(df.columns)}")
                        else:
                            # Try common column names
                            possible_columns = ['HGNC symbol', 'HGNC', 'Gene', 'Gene Symbol', 'Symbol', 'TF']
                            tf_column = None
                            for col_name in possible_columns:
                                if col_name in df.columns:
                                    tf_column = df[col_name]
                                    break
                            
                            if tf_column is None:
                                tf_column = df.iloc[:, 1] if len(df.columns) > 1 else df.iloc[:, 0]
                        
                        for tf in tf_column:
                            if pd.notna(tf):
                                tf_str = str(tf).strip().strip('"').strip("'").strip()
                                if tf_str:
                                    tfs.add(tf_str.upper())
                    else:
                        # Plain text (one TF per line)
                        for line in f:
                            tf = line.strip().strip('"').strip("'").strip()
                            if tf and not tf.startswith('#'):
                                tfs.add(tf.upper())
            
            elif file_ext == '.csv':
                df = pd.read_csv(tf_file, header=0 if header else None)
                if column:
                    tf_column = df[column]
                else:
                    tf_column = df.iloc[:, 0]
                
                for tf in tf_column:
                    if pd.notna(tf):
                        tf_str = str(tf).strip().strip('"').strip("'").strip()
                        if tf_str:
                            tfs.add(tf_str.upper())
            
            elif file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(tf_file, header=0 if header else None)
                if column:
                    tf_column = df[column]
                else:
                    tf_column = df.iloc[:, 0]
                
                for tf in tf_column:
                    if pd.notna(tf):
                        tf_str = str(tf).strip().strip('"').strip("'").strip()
                        if tf_str:
                            tfs.add(tf_str.upper())
            else:
                # Try as plain text
                with open(tf_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        tf = line.strip().strip('"').strip("'").strip()
                        if tf and not tf.startswith('#'):
                            tfs.add(tf.upper())
        
        except Exception as e:
            raise ValueError(f"Error reading TF file '{tf_file}': {e}")
        
        if self.verbose:
            print(f"  Loaded {len(tfs)} unique TFs")
        
        return tfs
    
    def load_motif_database(
        self,
        motif_file: str,
        tf_column: Optional[str] = None,
        motif_column: Optional[str] = None,
        header: bool = True
    ) -> Dict[str, Set[str]]:
        """
        Load motif database.
        
        Supports different formats:
        - CSV/TSV: Tabular format with TF and motif columns
        - Feather: pySCENIC cisTarget database format
        - Annotation files: Files mapping TFs to motifs
        
        Parameters:
        -----------
        motif_file : str
            Path to motif database file
        tf_column : str, optional
            Column name containing TF names
        motif_column : str, optional
            Column name containing motif names
        header : bool, default=True
            Whether file has a header row
        
        Returns:
        --------
        motif_dict : dict
            Dictionary mapping TF names (uppercase) to sets of motif names
        """
        if not os.path.exists(motif_file):
            raise FileNotFoundError(f"Motif file not found: {motif_file}")
        
        if self.verbose:
            print(f"Loading motif database from {motif_file}...")
        
        motif_dict: Dict[str, Set[str]] = {}
        file_ext = Path(motif_file).suffix.lower()
        
        try:
            if file_ext in ['.csv', '.tsv', '.txt']:
                # Tabular format
                sep = '\t' if file_ext == '.tsv' or (file_ext == '.txt' and '\t' in open(motif_file).readline()) else ','
                df = pd.read_csv(motif_file, sep=sep, header=0 if header else None)
                
                # Try to find TF and motif columns
                if tf_column and motif_column:
                    if tf_column in df.columns and motif_column in df.columns:
                        for _, row in df.iterrows():
                            tf = str(row[tf_column]).strip().upper()
                            motif = str(row[motif_column]).strip()
                            if pd.notna(tf) and pd.notna(motif) and tf and motif:
                                if tf not in motif_dict:
                                    motif_dict[tf] = set()
                                motif_dict[tf].add(motif)
                else:
                    tf_col = None
                    motif_col = None
                    
                    for col in df.columns:
                        col_str = str(col).upper()
                        if not tf_col and any(x in col_str for x in ['TF', 'TRANSCRIPTION', 'GENE', 'SYMBOL', 'HGNC']):
                            tf_col = col
                        if not motif_col and any(x in col_str for x in ['MOTIF', 'BINDING', 'TARGET']):
                            motif_col = col
                    
                    if tf_col and motif_col:
                        for _, row in df.iterrows():
                            tf = str(row[tf_col]).strip().upper()
                            motif = str(row[motif_col]).strip()
                            if pd.notna(tf) and pd.notna(motif) and tf and motif:
                                if tf not in motif_dict:
                                    motif_dict[tf] = set()
                                motif_dict[tf].add(motif)
                    else:
                        # Try first two columns
                        if len(df.columns) >= 2:
                            for _, row in df.iterrows():
                                tf = str(row.iloc[0]).strip().upper()
                                motif = str(row.iloc[1]).strip()
                                if pd.notna(tf) and pd.notna(motif) and tf and motif:
                                    if tf not in motif_dict:
                                        motif_dict[tf] = set()
                                    motif_dict[tf].add(motif)
            
            elif file_ext in ['.xlsx', '.xls']:
                # Excel format
                df = pd.read_excel(motif_file, header=0 if header else None)
                
                # Similar logic to CSV
                if tf_column and motif_column:
                    for _, row in df.iterrows():
                        tf = str(row[tf_column]).strip().upper()
                        motif = str(row[motif_column]).strip()
                        if pd.notna(tf) and pd.notna(motif) and tf and motif:
                            if tf not in motif_dict:
                                motif_dict[tf] = set()
                            motif_dict[tf].add(motif)
        
        except Exception as e:
            raise ValueError(f"Error reading motif file '{motif_file}': {e}")
        
        if self.verbose:
            print(f"  Loaded motifs for {len(motif_dict)} TFs")
            total_motifs = sum(len(motifs) for motifs in motif_dict.values())
            print(f"  Total TF-motif relationships: {total_motifs}")
        
        return motif_dict
    
    def calculate_overlap(
        self,
        tf_set: Set[str],
        motif_dict: Dict[str, Set[str]]
    ) -> Dict[str, any]:
        """
        Calculate overlap statistics between TF list and motif database.
        
        Parameters:
        -----------
        tf_set : set
            Set of TF gene symbols
        motif_dict : dict
            Dictionary mapping TF names to motif sets
        
        Returns:
        --------
        stats : dict
            Dictionary with overlap statistics
        """
        n_tfs = len(tf_set)
        tfs_in_database = tf_set.intersection(set(motif_dict.keys()))
        n_tfs_with_motifs = len(tfs_in_database)
        
        if n_tfs == 0:
            percent_with_motifs = 0.0
        else:
            percent_with_motifs = (n_tfs_with_motifs / n_tfs) * 100
        
        # Count total motifs
        total_motifs = sum(len(motifs) for tf, motifs in motif_dict.items() if tf in tf_set)
        avg_motifs_per_tf = total_motifs / n_tfs_with_motifs if n_tfs_with_motifs > 0 else 0
        
        missing_tfs = tf_set - set(motif_dict.keys())
        
        stats = {
            "n_tfs": n_tfs,  # Total number of TFs in the input TF list
            "n_tfs_with_motifs": n_tfs_with_motifs,  # Number of TFs from the list that have motifs in the database
            "percent_with_motifs": percent_with_motifs,  # Percentage of TFs in the list that have motifs (coverage percentage)
            "total_motifs": total_motifs,  # Total number of motif-TF relationships for TFs in the input list
            "avg_motifs_per_tf": avg_motifs_per_tf,  # Average number of motifs per TF (only counting TFs that have motifs)
            "missing_tfs": missing_tfs,  # Set of TFs from the input list that are NOT found in the motif database
            "tfs_with_motifs": tfs_in_database,  # Set of TFs from the input list that ARE found in the motif database
        }
        
        return stats
    
    def analyze_overlap(
        self,
        tf_file: str,
        motif_file: str,
        tf_column: Optional[str] = None,
        motif_tf_column: Optional[str] = None,
        motif_column: Optional[str] = None
    ) -> Dict[str, any]:
        """
        Analyze overlap between TF list and motif database.
        
        Parameters:
        -----------
        tf_file : str
            Path to TF list file
        motif_file : str
            Path to motif database file
        tf_column : str, optional
            Column name in TF file
        motif_tf_column : str, optional
            Column name for TFs in motif file
        motif_column : str, optional
            Column name for motifs in motif file
        
        Returns:
        --------
        results : dict
            Dictionary with overlap statistics
        """
        # Load TF list
        tf_set = self.load_tf_list(tf_file, column=tf_column)
        
        # Load motif database
        motif_dict = self.load_motif_database(
            motif_file,
            tf_column=motif_tf_column,
            motif_column=motif_column
        )
        
        # Calculate overlap
        stats = self.calculate_overlap(tf_set, motif_dict)
        
        return stats
    
    def print_summary(self, stats: Dict[str, any]):
        """
        Print a simple summary with just the percentage.
        
        Parameters:
        -----------
        stats : dict
            Statistics dictionary from calculate_overlap()
        """
        print(f"\nMotif Coverage: {stats['percent_with_motifs']:.2f}%")


def main():
    """Main function to run motif overlap analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Test overlap between TF list and motif database"
    )
    parser.add_argument(
        "--tf-file",
        type=str,
        help="Path to TF list file"
    )
    parser.add_argument(
        "--motif-file",
        type=str,
        help="Path to motif database file (.feather, .csv, .tsv, .xlsx)"
    )
    parser.add_argument(
        "--tf-column",
        type=str,
        default=None,
        help="Column name in TF file (for CSV/Excel)"
    )
    parser.add_argument(
        "--motif-tf-column",
        type=str,
        default=None,
        help="Column name for TFs in motif file"
    )
    parser.add_argument(
        "--motif-column",
        type=str,
        default=None,
        help="Column name for motifs in motif file"
    )
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = MotifOverlapAnalyzer(verbose=False)
    
    # Run analysis
    stats = analyzer.analyze_overlap(
        args.tf_file,
        args.motif_file,
        tf_column=args.tf_column,
        motif_tf_column=args.motif_tf_column,
        motif_column=args.motif_column
    )
    
    # Print summary
    analyzer.print_summary(stats)


if __name__ == "__main__":
    main()

