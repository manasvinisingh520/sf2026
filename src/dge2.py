"""
Differential Gene Expression (DGE) Analysis Class

This module provides a class-based interface for performing DGE analysis
using PyDESeq2 (DESeq2 implementation in Python).
"""

import pandas as pd
import numpy as np
from typing import Optional, Tuple, Dict, List, Any
from anndata import AnnData
from utils import filter_anndata_object

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
except ImportError:
    PYDESEQ2_AVAILABLE = False
    print("Warning: pydeseq2 not installed. Install with: pip install pydeseq2")


class DGEAnalyzer:
    
    def __init__(
        self,
        data_file: AnnData,
        condition: str,
        groups: List[Tuple[str, str]],
        method: str = "deseq2",
        threshold: Optional[Dict[str, float]] = None,
        topK: int = 10,
        output_file: str = "dge_results.csv",
        min_counts: int = 10,
        n_cpus: int = 4,
        control_genes: Optional[List[str]] = None
    ):
        """
        Initialize DGEAnalyzer.
        
        Parameters:
        -----------
        data_file : AnnData
            AnnData object with gene expression data (cells x genes format)
            Must have raw counts in .X and metadata in .obs
        condition : str
            Column name in adata.obs specifying the condition/group variable
        groups : list of tuples
            List of (group1, group2) tuples for multiple comparisons
            Example: [("treated", "control"), ("disease", "healthy")]
        method : str, optional
            DGE method to use (default: "deseq2")
        threshold : dict, optional
            Dictionary with 'padj' and 'log2fc' thresholds for significance
            Default: {'padj': 0.05, 'log2fc': 1.0}
        topK : int, optional
            Number of top genes to return (default: 10)
        output_file : str, optional
            Output CSV file name (default: "dge_results.csv")
        min_counts : int, optional
            Minimum total counts across samples to keep a gene (default: 10)
        n_cpus : int, optional
            Number of CPUs for parallel processing (default: 4)
        control_genes : list, optional
            List of control gene names for normalization (default: None)
        """
        # Validate pydeseq2 availability
        if method == "deseq2" and not PYDESEQ2_AVAILABLE:
            raise ImportError(
                "pydeseq2 is required for DESeq2 method. "
                "Install with: pip install pydeseq2"
            )
        
        # Validate condition column exists
        if condition not in data_file.obs.columns:
            raise ValueError(
                f"Condition column '{condition}' not found in adata.obs. "
                f"Available columns: {list(data_file.obs.columns)}"
            )
        
        # Validate and process groups
        if not groups:
            raise ValueError("'groups' list must contain at least one (group1, group2) tuple")
        
        unique_conditions = data_file.obs[condition].unique()
        
        # Validate all groups exist
        for g1, g2 in groups:
            if g1 not in unique_conditions:
                raise ValueError(
                    f"Group '{g1}' not found in condition '{condition}'. "
                    f"Available groups: {list(unique_conditions)}"
                )
            if g2 not in unique_conditions:
                raise ValueError(
                    f"Group '{g2}' not found in condition '{condition}'. "
                    f"Available groups: {list(unique_conditions)}"
                )
        
        # Store groups
        self.groups = groups        # Store parameters
        self.data_file = data_file
        self.condition = condition
        self.method = method
        self.topK = topK
        self.output_file = output_file
        self.min_counts = min_counts
        self.n_cpus = n_cpus
        self.control_genes = control_genes
        
        # Set default thresholds if not provided
        if threshold is None:
            threshold = {'padj': 0.05, 'log2fc': 1.0}
        self.threshold = threshold
        
        # Initialize results attributes
        self.dds = None  # DeseqDataSet object
        self.stat_res_dict = {}  # Dictionary of DeseqStats objects for each group pair
        self.results_dict = {}  # Dictionary of results DataFrames for each group pair
        self.significant_genes_dict = {}  # Dictionary of significant genes for each group pair
        
        # Store extracted counts and metadata
        self.counts_df = None
        self.metadata_df = None
    
    def _prepare_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Prepare count matrix and metadata from AnnData object.
        
        Returns:
        --------
        counts_df : pd.DataFrame
            Count matrix with samples (cells) as rows, genes as columns
            Index must match metadata_df index (cell names)
        metadata_df : pd.DataFrame
            Sample metadata with condition information
        """

        # Apply quality filters
        if self.min_genes is not None:
            self.data_file = filter_anndata_object(self.data_file, min_genes=self.min_genes, min_cells=3, min_counts=None, max_counts=None)
        
        if self.data_file is None:
            raise ValueError("Failed to filter AnnData object. Check your filter thresholds.")

        # Extract counts matrix (cells x genes from AnnData)
        # PyDESeq2's DeseqDataSet expects samples (cells) as rows, genes as columns
        # The row index (cell names) must match metadata_df index
        if hasattr(self.data_file.X, 'toarray'):
            # Sparse matrix
            counts_matrix = self.data_file.X.toarray()  # Keep as cells x genes
        else:
            # Dense matrix
            counts_matrix = self.data_file.X  # Keep as cells x genes
        
        # Create counts DataFrame with cell names as index (rows), gene names as columns
        self.counts_df = pd.DataFrame(
            counts_matrix,
            index=self.data_file.obs_names,  # Cells as rows - matches metadata_df index
            columns=self.data_file.var_names  # Genes as columns
        )
        
        # Extract metadata
        self.metadata_df = self.data_file.obs.copy()
        
        # Ensure indices match (important for PyDESeq2)
        # Filter metadata to only include cells present in counts
        # self.metadata_df = self.metadata_df.loc[self.counts_df.index] (CHECK)
        
        # Filter lowly expressed genes
        if 0:
            if self.min_counts > 0:
                gene_sums = self.counts_df.sum(axis=0)
                self.counts_df = self.counts_df.loc[:, gene_sums >= self.min_counts]
            
            # Verify indices match
            if not self.counts_df.index.equals(self.metadata_df.index):
                raise ValueError("Index mismatch between counts and metadata")
        
        return self.counts_df, self.metadata_df
    
    def run_dge(self) -> pd.DataFrame:
        """
        Run differential gene expression analysis.
        
        Returns:
        --------
        results : pd.DataFrame
            DataFrame with DGE results containing columns:
            - baseMean: Average normalized count
            - log2FoldChange: Log2 fold change (group1 vs group2)
            - lfcSE: Standard error of log2FoldChange
            - stat: Wald statistic
            - pvalue: Raw p-value
            - padj: Adjusted p-value (FDR)
        """
        if self.method != "deseq2":
            raise NotImplementedError(f"Method '{self.method}' not yet implemented. Only 'deseq2' is supported.")
        
        # Prepare data
        counts_df, metadata_df = self._prepare_data()
        
        # Create design formula
        design = f"~{self.condition}"
        
        # Initialize inference
        inference = DefaultInference()
        
        # Create DeseqDataSet
        self.dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design=design,
            inference=inference,
            n_cpus=self.n_cpus,
            control_genes=self.control_genes
        )
        
        # Run DESeq2 pipeline
        self.dds.deseq2()
        
        return self.dds
    
    def get_results(
        self,
        groups: Optional[List[Tuple[str, str]]] = None,
        save: bool = True,
        filename_prefix: Optional[str] = None
    ) -> Dict[str, Dict[str, Any]]:
        """
        Extract results, find significant genes, get top genes, and save results
        for all group pairs.
        
        This method combines:
        - Extracting results using DeseqStats (lines 251-258 equivalent)
        - Finding significant genes
        - Getting top genes
        - Saving results to files
        
        Parameters:
        -----------
        groups : list of tuples, optional
            List of (group1, group2) tuples for comparisons.
            If None, uses self.groups from initialization (default: None)
        save : bool, optional
            Whether to save results to files (default: True)
        filename_prefix : str, optional
            Prefix for output filenames. If None, uses self.output_file as base (default: None)
            
        Returns:
        --------
        results_dict : dict
            Dictionary with keys as f"{group1}_vs_{group2}" containing:
            - 'stat_res': DeseqStats object
            - 'results': All results DataFrame
            - 'significant_genes': Significant genes DataFrame
            - 'top_genes': Top K genes DataFrame
        """
        if self.dds is None:
            raise ValueError("DeseqDataSet not available. Run run_dge() first.")
        
        # Use provided groups or self.groups
        if groups is None:
            groups = self.groups
        
        if not groups:
            raise ValueError("No groups provided.")
        
        # Initialize inference
        inference = DefaultInference()
        
        # Dictionary to store all results
        all_results = {}
        
        # Iterate over all group pairs
        for group1, group2 in groups:
            stat_res = DeseqStats(
                self.dds,
                contrast=[self.condition, group1, group2],
                inference=inference,
                n_cpus=self.n_cpus
            )
            
            stat_res.summary()
            results = stat_res.results_df
            
            # Store stat_res and results
            pair_key = f"{group1}_vs_{group2}"
            self.stat_res_dict[pair_key] = stat_res
            self.results_dict[pair_key] = results
            
            # Find significant genes
            padj_threshold = self.threshold.get('padj', 0.05)
            log2fc_threshold = self.threshold.get('log2fc', 1.0)
            
            mask = (
                (results['padj'] < padj_threshold) &
                (abs(results['log2FoldChange']) > log2fc_threshold)
            )
            
            significant_genes = results[mask].copy()
            self.significant_genes_dict[pair_key] = significant_genes
            
            # Get genes with padj < 0.05
            padj_mask = (results['padj'] < 0.05) & (results['padj'].notna())
            genes_padj_005 = results[padj_mask].copy()
            if len(genes_padj_005) > 0:
                genes_padj_005 = genes_padj_005.sort_values('padj', ascending=True)
            
            # Keep top_genes
            top_genes = genes_padj_005.head(self.topK) if len(genes_padj_005) > 0 else pd.DataFrame()
            
            # Save results if requested
            if save:
                if filename_prefix is None:
                    base = self.output_file.rsplit('.', 1)[0] if '.' in self.output_file else self.output_file
                    if pair_key in base:
                        filename = f"{base}.csv"
                    else:
                        filename = f"{base}_{pair_key}.csv"
                else:
                    filename = f"{filename_prefix}_{pair_key}.csv"
                
                results.to_csv(filename)
            
            # Store all results for this pair
            all_results[pair_key] = {
                'stat_res': stat_res,
                'results': results,
                'significant_genes': significant_genes,
                'genes_padj_005': genes_padj_005,
                'top_genes': top_genes
            }
        
        return all_results
    
    def run_full_analysis(self) -> Dict[str, Any]:
        """
        Run complete DGE analysis pipeline.
        
        This method runs all steps:
        1. Run DGE analysis
        2. Get results for all group pairs (extracts, filters, gets top genes, saves)
        
        Returns:
        --------
        results_dict : dict
            Dictionary with keys as f"{group1}_vs_{group2}" containing:
            - 'stat_res': DeseqStats object
            - 'results': All DGE results
            - 'significant_genes': Significant DEGs
            - 'top_genes': Top K genes
        """
        # Step 1: Run DGE
        self.run_dge()
        
        # Step 2: Get results for all group pairs
        all_results = self.get_results()
        
        return all_results
    
    def __repr__(self):
        """String representation of the DGEAnalyzer object."""
        status = "Analysis run" if len(self.results_dict) > 0 else "Not run"
        if len(self.groups) > 1:
            comparisons = f"{len(self.groups)} comparisons: {self.groups}"
        else:
            comparisons = f"{self.group1} vs {self.group2}"
        return (
            f"DGEAnalyzer(\n"
            f"  method='{self.method}',\n"
            f"  condition='{self.condition}',\n"
            f"  comparisons={comparisons},\n"
            f"  threshold={self.threshold},\n"
            f"  topK={self.topK},\n"
            f"  output_file='{self.output_file}',\n"
            f"  status='{status}'\n"
            f")"
        )


def shuffle_data(adata: AnnData, seed: int) -> AnnData:
    rng = np.random.RandomState(seed)
    n_cells = adata.n_obs
    shuffled_indices = rng.permutation(n_cells)
    return adata[shuffled_indices].copy()


# Example usage
if __name__ == "__main__":
    """
    Example usage of DGEAnalyzer class.
    
    This example demonstrates how to use the class with an AnnData object
    loaded from MTX files (same as read_mtx_example.py).
    
    Usage:
        python dge.py --group1 GROUP1 --group2 GROUP2 --region REGION
        
    Example:
        python dge.py --group1 "Control" --group2 "Treatment" --region "EC"
    """
    import argparse
    import os
    
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Run Differential Gene Expression (DGE) analysis using PyDESeq2',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument(
        '--region',
        type=str,
        required=True,
        help='Region identifier to locate data files (e.g., EC, ITG, PFC, V1, V2)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=100,
        help='Random seed for cell shuffling before pseudobulk binning (default: 100)'
    )
    parser.add_argument(
        '--topK',
        type=int,
        default=10,
        help='Number of top genes to display (default: 50)'
    )
    parser.add_argument(
        '--bins',
        type=int,
        default=800,
        help='Target number of cells per bin for pseudobulk aggregation (default: 800)'
    )
    
    # Parse arguments
    args = parser.parse_args()
    random_seed = args.seed
    
    # Set global random seed for reproducibility
    np.random.seed(random_seed)
    
    # Import utilities and config
    from utils import read_mtx_file, create_anndata_object, read_excel_columns, get_region_file_paths
    from config import REGION_TO_TAB, DEFAULT_METADATA_PATH, DATA_DIR, METADATA_DATE_PREFIX
    
    # Configuration
    base_dir = DATA_DIR
    date_prefix = METADATA_DATE_PREFIX
    region = args.region
    
    # Get file paths using utility function
    mtx_path, row_annotation_path, col_annotation_path = get_region_file_paths(
        region,
        data_dir=base_dir,
        base_prefix=f"{date_prefix}_Astrocytes_{{region}}"
    )
    metadata_path = os.path.join(base_dir, f"{date_prefix}_Astrocytes_metadata.xlsx")
    
    # Get tab index from config
    tab_index = REGION_TO_TAB.get(region)
    if tab_index is None:
        raise ValueError(f"Region '{region}' not mapped to a tab index. Available regions: {list(REGION_TO_TAB.keys())}")
    
    # Load metadata
    metadata = read_excel_columns(metadata_path, columns=['cell_annotation', "RIN", "Path..Group.", "SampleName", 
                                                          "Median.UMI.Counts.per.Cell", "percent.mito"],
                                  sheet_name=tab_index)
    
    # Read the MTX file
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=str(mtx_path),
        row_annotation_path=str(row_annotation_path),
        col_annotation_path=str(col_annotation_path),
        transpose=False
    )
    
    # Create AnnData object
    adata = create_anndata_object(
        matrix=matrix,
        gene_names=gene_names,
        cell_names=cell_names,
        transpose=True,  # Matrix is genes x cells, transpose for AnnData (cells x genes)
        obs=metadata,
    )

    print (f"Before filtering: {adata.shape} (cells × genes)")

    # Apply quality filters
    adata = filter_anndata_object(adata, min_genes=200, min_cells=10, min_counts=None, max_counts=None, mito_max=0.15)
    
    print (f"After filtering: {adata.shape} (cells × genes)")

    if adata is None:
        raise ImportError("Failed to create AnnData object. Make sure anndata is installed.")
    
    # Determine condition column
    if 'Path..Group.' not in adata.obs.columns:
        raise ValueError(f"No condition column found. Available columns: {list(adata.obs.columns)}")
    
    condition_col = 'Path..Group.'
    
    # CRITICAL: Aggregate cells into pseudobulk samples if dataset is too large
    # DESeq2 is designed for bulk RNA-seq with few samples (3-10 per condition)
    # Single-cell data with thousands of cells causes memory errors
    # Pseudobulk aggregation is standard practice for single-cell DGE analysis
    max_cells_for_deseq2 = 500  # DESeq2 can't handle more than ~500 cells efficiently
    target_cells_per_bin = args.bins  # Target number of cells per bin
    
    if n_after > max_cells_for_deseq2:
        from scipy.sparse import csr_matrix
        
        pseudobulk_data = []
        pseudobulk_metadata = []
        
        for group_name in [group1, group2]:
            group_mask = adata_filtered.obs[condition_col] == group_name
            group_data = adata_filtered[group_mask]
            unique_patients = group_data.obs['SampleName'].dropna().unique()
            
            print(f"\n{group_name}:")
            total_bins_group = 0
            
            for patient_name in unique_patients:
                patient_mask = group_data.obs['SampleName'] == patient_name
                patient_cells = group_data[patient_mask]
                n_cells_per_patient = patient_cells.n_obs
                
                if n_cells_per_patient < 80:
                    print(f"  {patient_name}: {n_cells_per_patient} cells (skipped, < 80 cells)")
                    continue
                
                # Shuffle cells for this patient using the seed from parser argument
                if args.seed != 100:
                    patient_cells = shuffle_data(patient_cells, random_seed)
                
                n_bins_for_patient = max(1, int(np.ceil(n_cells_per_patient / target_cells_per_bin)))
                total_bins_group += n_bins_for_patient
                print(f"  {patient_name}: {n_cells_per_patient} cells {n_bins_for_patient} bins")
                base_cells_per_bin = n_cells_per_patient // n_bins_for_patient
                remainder = n_cells_per_patient % n_bins_for_patient
                
                cell_idx = 0
                for bin_idx in range(n_bins_for_patient):
                    bin_size = base_cells_per_bin + (1 if bin_idx < remainder else 0)
                    start_idx = cell_idx
                    end_idx = cell_idx + bin_size
                    
                    patient_cells_subset = patient_cells[start_idx:end_idx]
                    cell_idx = end_idx
                    
                    if hasattr(patient_cells_subset.X, 'toarray'):
                        bin_counts = patient_cells_subset.X.sum(axis=0).A1
                    else:
                        bin_counts = patient_cells_subset.X.sum(axis=0)
                    
                    pseudobulk_data.append(bin_counts)
                    pseudobulk_metadata.append({
                        condition_col: group_name,
                        'SampleName': patient_name,
                        '_sample_id': f"{group_name}_{patient_name}_bin{bin_idx+1}"
                    })
            
            print(f"  Total bins for {group_name}: {total_bins_group}")
        
        print(f"\nTotal pseudobulk samples: {len(pseudobulk_data)}")
        pseudobulk_matrix = np.array(pseudobulk_data)
        pseudobulk_obs = pd.DataFrame(pseudobulk_metadata)
        pseudobulk_obs.index = pseudobulk_obs['_sample_id']
        
        adata_filtered = AnnData(
            X=csr_matrix(pseudobulk_matrix),
            obs=pseudobulk_obs,
            var=adata_filtered.var.copy()
        )
    
    # Create output directory if it doesn't exist
    output_dir = "results/dge4"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create DGE analyzer with group comparisons (using filtered/aggregated data)
    output_file = os.path.join(output_dir, f"dge_results_{region}_{group1}_vs_{group2}_bins{target_cells_per_bin}_seed{random_seed}.csv")
    dge = DGEAnalyzer(
        data_file=adata_filtered,
        condition=condition_col,
        groups=[(group1, group2)],
        method="deseq2",
        threshold={'padj': 0.05, 'log2fc': 0.5},
        topK=args.topK,
        output_file=output_file,
        min_counts=10,
        n_cpus=2
    )
    
    # Run full analysis
    results = dge.run_full_analysis()
    
    # Access results for each comparison
    print("\n" + "=" * 60)
    print("Results Summary")
    print("=" * 60)
    for comparison_key, comparison_results in results.items():
        print(f"\nComparison: {comparison_key}")
        print(f"Total genes tested: {len(comparison_results['results'])}")
        print(f"Significant DEGs (padj < 0.05, |log2FC| > 0.5): {len(comparison_results['significant_genes'])}")
        
        genes_padj_005 = comparison_results.get('genes_padj_005', pd.DataFrame())
        top_genes = comparison_results.get('top_genes', pd.DataFrame())
        
        print(f"\nGenes with padj < 0.05: {len(genes_padj_005)}")
        if len(top_genes) > 0:
            print(f"\nTop {args.topK} genes (by padj):")
            print(top_genes[['log2FoldChange', 'padj', 'baseMean']])
        elif len(genes_padj_005) > 0:
            # Fallback: show top K from genes_padj_005 if top_genes is empty
            print(f"\nTop {args.topK} genes (by padj):")
            print(genes_padj_005.head(args.topK)[['log2FoldChange', 'padj', 'baseMean']])
        else:
            print("  No genes found with padj < 0.05")

