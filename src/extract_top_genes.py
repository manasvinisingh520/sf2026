"""
Extract genes with padj < 0.05 from DGE results files.

This script reads DGE results CSV files and prints all genes with padj < 0.05.
"""

import pandas as pd
import argparse
import os
from pathlib import Path
from utils import get_dge_results_dir
from config import REGIONS

# Hardcoded gene list - add your genes here
GENE_LIST = [
    "ANGPTL4", "FGF2", "HGF", "NRP1", "NTRK3", "TGFB2", "TGFBR3", "VEGFA",
    "LGALS3", "SCARA3", "MAP3K14", "MAP4K4", "GPC4", "MATN2", "SERPINA3",
    "SERPINE2", "RFTN1",

    "DPY19L3", "EXT1", "GALNT2", "HIF1A", "MAN2A1", "ST6GALNAC3", "SULF1",
    "UGP2", "CHI3L1", "IFI16", "IL1R1", "IL6R", "PLA2G4C",

    "ACSL3", "APPL2", "ELOVL5", "ELOVL6", "GPCPD1", "LDLR", "LPIN1",
    "NCEH1", "PLPP3", "SGMS2", "SLC44A3", "PPP5C", "SPRED1", "SPRED2",

    "BAG3", "CRYAB", "HSP90AA1", "HSPA1A", "HSPB1", "HSPH1",
    "ATP5F1E", "ATP5ME", "COX6A1", "COX6C", "NDUFB2", "NDUFA4", "UQCRB",
    "MT1E", "MT1F", "MT1G", "MT1M", "MT1X", "MT3",
    "NFE2L2", "PRDX1", "SOD1", "SOD2", "IL17RB", "IRS2", "STAT3",
    "IFITM2", "IFITM3",
    "EEF1A1", "EIF2S2", "RPL37A", "RPL38", "RPLP2", "RPS21", "RPS24", "RPS28",

    "CD47", "DPP10", "GABRB1", "GRIA2", "ITPR2",
    "AASS", "ALDH1A1", "GLUD1", "MTR", "SLC1A3",
    "KCNIP1", "KCNQ3", "NALCN", "RYR3", "SCN11A",
    "SHISA6", "SLC13A3", "SLC4A4",
    "ATP9B", "TRPC1", "UNC79", "UNC80",

    "CRB1", "GFAP", "GLUL", "SLC1A2",
    "ACACB", "ACSS1", "CHKA", "DGKG", "DKK3", "MALRD1",
    "PLCE1", "SLC44A1", "SREBF1",
    "CDH20", "FRMPD2", "ITGB4", "NLGN4Y", "NTRK2", "SDC4", "TESK2",
    "ATP13A4", "CACNA2D3", "CNTN1", "KCNN3", "NKAIN2",
    "SLC24A4", "TPCN1", "TRPM3",

    "DNAJB2", "EEF2K", "HSPB8", "WFS1",
    "CPEB3", "CTSD", "KAT2A", "LRP4", "NEDD4L",
    "ATP2B4", "CALM1", "CKB", "GJA1", "KCNJ10", "SLC39A12",
    "COL27A1", "COL5A3", "PLEC",
    "APOE", "RORA", "RXRA", "SLC27A1"
]
"""

GENE_LIST = [
    "APOE", "APP", "PSEN1", "PSEN2", "ABCA7", "PICALM", "PLD3", "TREM2", "SORL1", "INPP5D", "CASS4", "NME8", "MEF2C",
    "PTK2B", "FERMT2", "ZCWPW1", "DSG2", "UNC5C", "ADAM10", "SLC24A4", "RIN3", "HLA-DRB1", "HLA-DRB5", "CELF1", 
    "PLD3", "AKAP9", "ATXN1", "CD33", "GWA14Q34", "DLGAP1", "CLU", "CR1", "BIN1", "CD2AP", "MS4A", "EPHA1"
]
"""


def compare_to_gene_list(genes_df: pd.DataFrame, label: str = "", return_normalized_names: bool = False, verbose: bool = True) -> tuple:
    """
    Compare genes in DataFrame to GENE_LIST.
    
    Parameters
    ----------
    genes_df : pd.DataFrame
        DataFrame with genes as index (gene names)
    label : str, optional
        Label for display purposes (e.g., "top 10", "padj < 0.05 results")
    print_results : bool, optional
        Whether to print comparison results (default: True)
    return_normalized_names : bool, optional
        Whether to return normalized overlapping gene names (default: False)
    
    Returns
    -------
    tuple: (num_matching, overlap_percentage, overlapping_genes_normalized)
        - num_matching: count of matching genes
        - overlap_percentage: percentage of GENE_LIST that was found
        - overlapping_genes_normalized: list of normalized matching gene names (if return_normalized_names=True)
    """
    # Normalize gene names for case-insensitive comparison
    normalize = lambda g: str(g).strip().upper()
    gene_list_normalized = [normalize(g) for g in GENE_LIST]
    genes_from_df_normalized = {normalize(g) for g in genes_df.index}
    
    # Find overlapping genes
    overlapping_genes_normalized = [g for g in gene_list_normalized if g in genes_from_df_normalized]
    num_matching = len(overlapping_genes_normalized)
    overlap_pct = (num_matching / len(GENE_LIST)) * 100
    
    if verbose:
        title = "Gene List Overlap Analysis"
        print(f"\n{'=' * 60}\n{title}\n{'=' * 60}")
        print(f"Overlap percentage: {overlap_pct:.2f}% ({num_matching}/{len(GENE_LIST)})")
        if num_matching > 0:
            display_cols = [col for col in ['log2FoldChange', 'padj', 'baseMean'] if col in genes_df.columns]
            print(f"\nMatching genes:")
            for norm_gene in overlapping_genes_normalized:
                print(f"{norm_gene}")
                print(genes_df.loc[norm_gene][display_cols])
                print()
    
    return (num_matching, overlap_pct, overlapping_genes_normalized) if return_normalized_names else (num_matching, overlap_pct)

def get_top_genes(results: pd.DataFrame, field_name: str, n: int = 10, padj_threshold: float = None, log2fc_threshold: float = None) -> pd.DataFrame:
    """
    Get the top n genes by a given field from a results DataFrame.

    Parameters
    ----------
    results : pd.DataFrame
        DataFrame with genes as index (gene names)
    field_name : str
        Field name to sort by
    n : int, optional
        Number of top genes to return (default: 10)
    """
    display_cols = [col for col in ['log2FoldChange', 'padj', 'baseMean'] if col in results.columns]
    if padj_threshold is not None:
        results = results[results['padj'] < padj_threshold]
    if log2fc_threshold is not None:
        results = results[abs(results['log2FoldChange']) > log2fc_threshold]
    
    short_type = True if field_name == 'padj' else False
    x = results.sort_values(field_name, ascending=short_type).head(min(n, len(results)))
    print (f"Top {len(x)} genes sorted by {field_name} with padj_max: {padj_threshold} and log2fc_min: {log2fc_threshold}:")
    print(x[display_cols])
    print()
    return x

def extract_top_genes(results_file: str = None, output_file: str = None, padj_threshold: float = 0.05, log2fc_threshold: float = 0.5, args: argparse.Namespace = None, combined_results: pd.DataFrame = None) -> pd.DataFrame:
    """
    Extract all genes with padj < threshold from DGE results file.
    Always prints top 10 genes passing padj threshold and top 10 passing log2fc threshold.
    
    Parameters
    ----------
    results_file : str
        Path to DGE results CSV file
    output_file : str, optional
        Path to save results. If None, only prints to stdout
    padj_threshold : float, optional
        padj threshold for significant genes (default: 0.05)
    log2fc_threshold : float, optional
        log2fc threshold for significant genes (default: 0.5)
    """
    if combined_results is not None:
        results = combined_results
        print(f"Using combined results from union of files")
    else:
        print(f"Reading DGE results from: {results_file}")
        results = pd.read_csv(results_file, index_col=0)
    print(f"Total genes in results: {len(results)}")
    
    # Print top based on padj and log2fc thresholds
    _ = get_top_genes(results, 'padj', args.top_n)
    _ = get_top_genes(results, 'log2FoldChange', args.top_n)
    _ = get_top_genes(results, 'padj', args.top_n, padj_threshold=args.padj_max)
    _ = get_top_genes(results, 'log2FoldChange', args.top_n, log2fc_threshold=args.log2fc_min)
    x = get_top_genes(results, 'padj', args.top_n, padj_threshold=args.padj_max, log2fc_threshold=args.log2fc_min)
    
    # Compare overlap with GENE_LIST for top 10, top 200, and top 500
    print(f"\n{'=' * 60}")
    print("GENE_LIST Overlap Analysis for Different Top N Values")
    print(f"{'=' * 60}")
    
    # Filter results by thresholds first
    filtered_results = results.copy()
    if args.padj_max is not None:
        filtered_results = filtered_results[filtered_results['padj'] < args.padj_max]
    if args.log2fc_min is not None:
        filtered_results = filtered_results[abs(filtered_results['log2FoldChange']) > args.log2fc_min]
    
    # Sort by padj and get top N for each
    filtered_results = filtered_results.sort_values('padj', ascending=True)
    
    for top_n in [10, 200, 500]:
        top_genes = filtered_results.head(min(top_n, len(filtered_results)))
        num_matching, overlap_pct, overlapping_genes = compare_to_gene_list(
            top_genes, f'top {top_n}', return_normalized_names=True, verbose=False
        )
        print(f"Top {top_n}: {overlap_pct:.2f}% overlap ({num_matching}/{len(GENE_LIST)} genes)", end="")
        if overlapping_genes:
            print(f" - Matching: {', '.join(overlapping_genes)}")
        else:
            print()
    
    return x


def main():
    parser = argparse.ArgumentParser(
        description='Extract all genes with padj < 0.05 from DGE results file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('-version', type=str, choices=['1', '2', '3', '4'], required=True, help='DGE version (1, 2, 3, or 4).')
    parser.add_argument('-region', type=str, choices=REGIONS, required=True, help='Region name (e.g., EC, ITG)')
    parser.add_argument('-cond1', type=str, required=True, help='Condition 1 identifier')
    parser.add_argument('-cond2', type=str, required=True, help='Condition 2 identifier')
    parser.add_argument('-bins', type=int, required=True, help='Bin number')
    parser.add_argument('-top_n', type=int, required=True, help='Number of top genes to extract. Uses same filtering as dge.py (padj < 0.05, |log2FC| > threshold)')
    parser.add_argument('-output_file', type=str, default=None, help='Output CSV file path to save results with padj < 0.05')
    parser.add_argument("-padj_max", type=float, default=0.05, help="padj threshold for significant genes")
    parser.add_argument("-log2fc_min", type=float, default=0.5, help="log2fc threshold for significant genes")
    args = parser.parse_args()
    
    # Get results directory using utility function
    version = int(args.version)
    base_dir = get_dge_results_dir(version)
    
    # Find all files matching region and bin size
    pattern = f"dge_results_{args.region}_*_bins{args.bins}*.csv"
    matching_files = sorted(list(base_dir.glob(pattern)))
    
    if len(matching_files) == 0:
        raise FileNotFoundError(
            f"No files found matching pattern: {pattern}\n"
            f"In directory: {base_dir.absolute()}"
        )
    
    print(f"Using DGE version: {version}")
    print(f"Searching in directory: {base_dir}")
    print(f"Found {len(matching_files)} files matching pattern: {pattern}")
    for f in matching_files:
        print(f"  - {f.name}")
    
    # Load all files and take union
    all_results = []
    for file_path in matching_files:
        df = pd.read_csv(file_path, index_col=0)
        all_results.append(df)
    
    # Combine all results - take union of genes, keep best padj for each gene
    combined_results = pd.concat(all_results, axis=0)
    # For genes that appear in multiple files, keep the one with lowest padj
    combined_results = combined_results.sort_values('padj', ascending=True, na_position='last')
    combined_results = combined_results[~combined_results.index.duplicated(keep='first')]
    
    print(f"\n{'=' * 60}")
    print(f"Union of {len(matching_files)} files: {len(combined_results)} unique genes")
    print(f"{'=' * 60}\n")
    
    # Extract genes with padj < threshold from union
    genes = extract_top_genes(
        results_file=None,  # We'll pass the combined results directly
        output_file=args.output_file,
        padj_threshold=args.padj_max,
        log2fc_threshold=args.log2fc_min,
        args=args,
        combined_results=combined_results
    )
    return genes


if __name__ == "__main__":
    main()

