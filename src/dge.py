"""
Differential Gene Expression (DGE) Analysis Class

This module provides a class-based interface for performing DGE analysis
using PyDESeq2 (DESeq2 implementation in Python).
"""

import pandas as pd
import numpy as np
from typing import Optional, Tuple, Dict, List, Any
from anndata import AnnData

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
except ImportError:
    PYDESEQ2_AVAILABLE = False
    print("Warning: pydeseq2 not installed. Install with: pip install pydeseq2")


class DGEAnalyzer:
    """
    Class for performing Differential Gene Expression (DGE) analysis.
    
    This class encapsulates the DGE workflow, providing methods to:
    - Run DGE analysis using DESeq2
    - Find significant differentially expressed genes
    - Extract top K genes by significance
    - Save results to file

    Methods:
    - run_dge(): Run DGE analysis
    - get_results(): Extract results for all group pairs
    """
    
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
            AnnData object with gene expression data (cells × genes format)
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
        print("Preparing data from AnnData object...")
        
        # Extract counts matrix (cells × genes from AnnData)
        # PyDESeq2's DeseqDataSet expects samples (cells) as rows, genes as columns
        # The row index (cell names) must match metadata_df index
        if hasattr(self.data_file.X, 'toarray'):
            # Sparse matrix
            counts_matrix = self.data_file.X.toarray()  # Keep as cells × genes
        else:
            # Dense matrix
            counts_matrix = self.data_file.X  # Keep as cells × genes
        
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
        self.metadata_df = self.metadata_df.loc[self.counts_df.index]
        
        # Filter lowly expressed genes (sum across samples/cells, axis=0)
        if self.min_counts > 0:
            gene_sums = self.counts_df.sum(axis=0)  # Sum across cells (rows)
            n_genes_before = len(self.counts_df.columns)
            self.counts_df = self.counts_df.loc[:, gene_sums >= self.min_counts]
            n_genes_after = len(self.counts_df.columns)
            print(f"Filtered genes: {n_genes_before} -> {n_genes_after} "
                  f"(kept genes with >= {self.min_counts} total counts)")
        
        print(f"Counts matrix shape: {self.counts_df.shape} (cells × genes)")
        print(f"Metadata shape: {self.metadata_df.shape}")
        print(f"Condition '{self.condition}': {self.metadata_df[self.condition].value_counts().to_dict()}")
        
        # Verify indices match
        if not self.counts_df.index.equals(self.metadata_df.index):
            raise ValueError(
                f"Index mismatch: counts_df index ({len(self.counts_df.index)} cells) "
                f"does not match metadata_df index ({len(self.metadata_df.index)} cells)"
            )
        
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
        print(f"\nDesign formula: {design}")
        
        # Initialize inference
        inference = DefaultInference()
        
        # Create DeseqDataSet
        print("Creating DeseqDataSet...")
        self.dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design=design,
            inference=inference,
            n_cpus=self.n_cpus,
            control_genes=self.control_genes
        )
        
        # Run DESeq2 pipeline
        print("Running DESeq2 pipeline...")
        self.dds.deseq2()
        print("DESeq2 pipeline complete.")
        
        print(f"\nDGE analysis complete. Use get_results() to extract results for group pairs.")
        
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
            raise ValueError("No groups provided. Specify groups in __init__ or get_results().")
        
        # Initialize inference
        inference = DefaultInference()
        
        # Dictionary to store all results
        all_results = {}
        
        print("=" * 60)
        print(f"Processing {len(groups)} group comparison(s)")
        print("=" * 60)
        
        # Iterate over all group pairs
        for idx, (group1, group2) in enumerate(groups, 1):
            print(f"\n[{idx}/{len(groups)}] Processing: {group1} vs {group2}")
            print("-" * 60)
            
            # Extract results (equivalent to lines 251-258)
            print(f"Extracting results: {group1} vs {group2}")
            stat_res = DeseqStats(
                self.dds,
                contrast=[self.condition, group1, group2],
                inference=inference,
                n_cpus=self.n_cpus
            )
            
            # Generate summary
            stat_res.summary()
            
            # Get results DataFrame
            results = stat_res.results_df
            print(f"Total genes tested: {len(results)}")
            
            # Store stat_res and results
            pair_key = f"{group1}_vs_{group2}"
            self.stat_res_dict[pair_key] = stat_res
            self.results_dict[pair_key] = results
            
            # Find significant genes
            padj_threshold = self.threshold.get('padj', 0.05)
            log2fc_threshold = self.threshold.get('log2fc', 1.0)
            
            print(f"\nFiltering significant genes...")
            print(f"  Thresholds: padj < {padj_threshold}, |log2FC| > {log2fc_threshold}")
            
            mask = (
                (results['padj'] < padj_threshold) &
                (abs(results['log2FoldChange']) > log2fc_threshold)
            )
            
            significant_genes = results[mask].copy()
            self.significant_genes_dict[pair_key] = significant_genes
            
            print(f"  Significant DEGs found: {len(significant_genes)}")
            
            # Print summary by direction
            if len(significant_genes) > 0:
                upregulated = (significant_genes['log2FoldChange'] > 0).sum()
                downregulated = (significant_genes['log2FoldChange'] < 0).sum()
                print(f"    Upregulated ({group1} > {group2}): {upregulated}")
                print(f"    Downregulated ({group1} < {group2}): {downregulated}")
            
            # Get top genes
            genes_df = significant_genes if len(significant_genes) > 0 else results
            sort_by = 'padj'
            ascending = True if sort_by in ['padj', 'pvalue'] else False
            top_genes = genes_df.sort_values(sort_by, ascending=ascending).head(self.topK)
            
            print(f"\nTop {self.topK} genes sorted by {sort_by}:")
            print(top_genes[['log2FoldChange', 'padj', 'baseMean']].head(self.topK))
            
            # Save results if requested
            if save:
                if filename_prefix is None:
                    # Extract base filename without extension
                    base = self.output_file.rsplit('.', 1)[0] if '.' in self.output_file else self.output_file
                    filename = f"{base}_{pair_key}.csv"
                else:
                    filename = f"{filename_prefix}_{pair_key}.csv"
                
                results.to_csv(filename)
                print(f"\nResults saved to '{filename}'")
            
            # Store all results for this pair
            all_results[pair_key] = {
                'stat_res': stat_res,
                'results': results,
                'significant_genes': significant_genes,
                'top_genes': top_genes
            }
        
        print("\n" + "=" * 60)
        print("All comparisons complete!")
        print("=" * 60)
        
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
        print("=" * 60)
        print("Running Full DGE Analysis Pipeline")
        print("=" * 60)
        
        # Step 1: Run DGE
        self.run_dge()
        
        # Step 2: Get results for all group pairs
        # This combines: extracting results, finding significant genes,
        # getting top genes, and saving results
        all_results = self.get_results()
        
        print("\n" + "=" * 60)
        print("Analysis Complete!")
        print("=" * 60)
        
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
        epilog="""
Examples:
  # Basic usage with two groups and region
  python dge.py --group1 "Control" --group2 "Treatment" --region "EC"
  
  # Different regions: EC, ITG, PFC, V1, V2
  python dge.py --group1 "A" --group2 "B" --region "V1"
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--group1',
        type=str,
        required=True,
        help='First group name for comparison (numerator in fold change)'
    )
    parser.add_argument(
        '--group2',
        type=str,
        required=True,
        help='Second group name for comparison (denominator in fold change)'
    )
    parser.add_argument(
        '--region',
        type=str,
        required=True,
        help='Region identifier to locate data files (e.g., EC, ITG, PFC, V1, V2)'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Import utilities for reading MTX files
    from utils import read_mtx_file, create_anndata_object, read_excel_columns
    
    # Construct file paths based on region
    base_dir = r"i:\sf2026\data"
    date_prefix = "2025-11-16"  # Adjust if needed
    region = args.region
    
    # Construct paths for region-specific files
    mtx_path = os.path.join(base_dir, f"{date_prefix}_Astrocytes_{region}_matrix.mtx")
    row_annotation_path = os.path.join(base_dir, f"{date_prefix}_Astrocytes_{region}_row_annotation.txt")
    col_annotation_path = os.path.join(base_dir, f"{date_prefix}_Astrocytes_{region}_cell_annotation.txt")
    metadata_path = os.path.join(base_dir, f"{date_prefix}_Astrocytes_metadata.xlsx")
    
    print("=" * 60)
    print("Loading Data from MTX Files")
    print("=" * 60)
    print(f"Region: {region}")
    print(f"Group 1: {args.group1}")
    print(f"Group 2: {args.group2}")
    print()
    
    # Load metadata
    metadata = read_excel_columns(metadata_path, columns=['cell_annotation', "RIN", "Path..Group."])
    print(f"Loaded metadata: {metadata.shape}")
    
    # Read the MTX file
    matrix, gene_names, cell_names = read_mtx_file(
        mtx_path=mtx_path,
        row_annotation_path=row_annotation_path,
        col_annotation_path=col_annotation_path,
        transpose=False  # Matrix will be genes × cells
    )
    
    print(f"Loaded matrix shape: {matrix.shape} (genes × cells)")
    
    # Filter metadata to include only rows where 'cell_annotation' matches the loaded cell_names
    if cell_names and 'cell_annotation' in metadata.columns:
        filtered_metadata = metadata[metadata['cell_annotation'].isin(cell_names)].copy()
        print(f"Filtered metadata shape: {filtered_metadata.shape}")
        # Align metadata index with cell_names for proper AnnData integration
        # Set cell_annotation as index and reindex to match cell_names order
        filtered_metadata = filtered_metadata.set_index('cell_annotation').reindex(cell_names)
        print(f"Cells with metadata: {filtered_metadata.index.notna().sum()}/{len(cell_names)}")
    else:
        filtered_metadata = metadata
        print("Note: 'cell_annotation' column not found in metadata or cell_names unavailable.")
    # Example: Access specific gene or cell
    print(f"Number of rows in filtered metadata: {len(filtered_metadata.index)}")
    
    # Create AnnData object
    adata = create_anndata_object(
        matrix=matrix,
        gene_names=gene_names,
        cell_names=cell_names,
        transpose=True,  # Matrix is genes × cells, transpose for AnnData (cells × genes)
        obs=filtered_metadata,
    )
    
    if adata is None:
        raise ImportError("Failed to create AnnData object. Make sure anndata is installed.")
    
    print(f"\nAnnData object created: {adata.shape} (cells × genes)")
    print(f"Available metadata columns: {list(adata.obs.columns)}")
    
    # Determine condition column (auto-detect Path..Group. or condition)
    if 'Path..Group.' in adata.obs.columns:
        condition_col = 'Path..Group.'
    else:
        raise ValueError(
            f"No condition column found. Available columns: {list(adata.obs.columns)}"
        )
    
    print(f"\nUsing '{condition_col}' as condition column")
    unique_conditions = adata.obs[condition_col].dropna().unique()
    print(f"Unique conditions in data: {list(unique_conditions)}")
    print(f"Condition counts:\n{adata.obs[condition_col].value_counts()}")
    
    # Use parsed groups
    group1 = args.group1
    group2 = args.group2
    
    # Validate that specified groups exist in the data
    if group1 not in unique_conditions:
        raise ValueError(
            f"Group '{group1}' not found in condition column '{condition_col}'. "
        )
    if group2 not in unique_conditions:
        raise ValueError(
            f"Group '{group2}' not found in condition column '{condition_col}'. "
        )
    
    # Check group sizes
    n_group1 = (adata.obs[condition_col] == group1).sum()
    n_group2 = (adata.obs[condition_col] == group2).sum()
    print(f"\nGroup sizes:")
    print(f"  {group1}: {n_group1} cells")
    print(f"  {group2}: {n_group2} cells")
    
    if n_group1 < 2 or n_group2 < 2:
        raise ValueError(
            f"Groups must have at least 2 samples each. "
            f"Found {group1}: {n_group1}, {group2}: {n_group2}"
        )
    
    # log2FC > 0 means higher in {group1}, log2FC < 0 means higher in {group2}
    print(f"\nComparing: {group1} vs {group2}")
    
    # CRITICAL: Filter AnnData to only include cells from the two groups being compared
    # This reduces memory usage significantly for large datasets
    print(f"\nFiltering AnnData to only include groups {group1} and {group2}...")
    mask = adata.obs[condition_col].isin([group1, group2])
    n_before = adata.n_obs
    adata_filtered = adata[mask].copy()
    n_after = adata_filtered.n_obs
    print(f"Filtered from {n_before:,} to {n_after:,} cells ({n_after/n_before*100:.1f}% of original)")
    
    # Check final group sizes
    group_counts_filtered = adata_filtered.obs[condition_col].value_counts()
    print(f"\nFinal group sizes after filtering:")
    print(f"  {group1}: {group_counts_filtered.get(group1, 0):,} cells")
    print(f"  {group2}: {group_counts_filtered.get(group2, 0):,} cells")
    
    # CRITICAL: Aggregate cells into pseudobulk samples if dataset is too large
    # DESeq2 is designed for bulk RNA-seq with few samples (3-10 per condition)
    # Single-cell data with thousands of cells causes memory errors
    # Pseudobulk aggregation is standard practice for single-cell DGE analysis
    max_cells_for_deseq2 = 500  # DESeq2 can't handle more than ~500 cells efficiently
    
    if True and n_after > max_cells_for_deseq2:
        print(f"\n⚠️  WARNING: Dataset has {n_after:,} cells - too large for DESeq2!")
        print(f"   DESeq2 will fail with memory errors above ~500 cells.")
        print(f"   Aggregating cells into pseudobulk samples...")
        
        from scipy.sparse import csr_matrix
        
        # Target: ~5-10 pseudobulk samples per group (standard for bulk RNA-seq)
        n_pseudobulk_per_group = 500
        pseudobulk_data = []
        pseudobulk_metadata = []
        
        for group_name in [group1, group2]:
            group_mask = adata_filtered.obs[condition_col] == group_name
            group_data = adata_filtered[group_mask]
            n_cells_in_group = group_data.n_obs
            
            # Calculate how many cells per pseudobulk sample
            cells_per_pseudobulk = max(1, n_cells_in_group // n_pseudobulk_per_group)
            
            # Create pseudobulk samples
            for i in range(n_pseudobulk_per_group):
                start_idx = i * cells_per_pseudobulk
                end_idx = min((i + 1) * cells_per_pseudobulk, n_cells_in_group)
                
                if start_idx >= n_cells_in_group:
                    break
                    
                sample_cells = group_data[start_idx:end_idx]
                
                # Sum counts across cells in this pseudobulk sample
                if hasattr(sample_cells.X, 'toarray'):
                    sample_counts = sample_cells.X.sum(axis=0).A1  # Sum across cells (axis=0)
                else:
                    sample_counts = sample_cells.X.sum(axis=0)
                
                pseudobulk_data.append(sample_counts)
                pseudobulk_metadata.append({
                    condition_col: group_name,
                    '_sample_id': f"{group_name}_pb{i+1}"
                })
        
        # Create new AnnData with pseudobulk samples
        pseudobulk_matrix = np.array(pseudobulk_data)
        pseudobulk_obs = pd.DataFrame(pseudobulk_metadata)
        pseudobulk_obs.index = pseudobulk_obs['_sample_id']
        
        adata_filtered = AnnData(
            X=csr_matrix(pseudobulk_matrix),
            obs=pseudobulk_obs,
            var=adata_filtered.var.copy()
        )
        
        print(f"✅ Successfully aggregated to {adata_filtered.n_obs} pseudobulk samples")
        print(f"   Sample distribution:")
        print(adata_filtered.obs[condition_col].value_counts())
        print(f"   Memory usage reduced by ~{n_after // adata_filtered.n_obs}x")
    else:
        print(f"\n✓ Dataset size ({n_after:,} cells) is acceptable for DESeq2.")
    
    # Create DGE analyzer with group comparisons (using filtered/aggregated data)
    dge = DGEAnalyzer(
        data_file=adata_filtered,
        condition=condition_col,
        groups=[(group1, group2)],
        method="deseq2",
        threshold={'padj': 0.05, 'log2fc': 0.5},
        topK=10,
        output_file=f"dge_results_{region}_{group1}_vs_{group2}.csv",
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
        print(f"Significant DEGs: {len(comparison_results['significant_genes'])}")
        print(f"\nTop {dge.topK} genes:")
        print(comparison_results['top_genes'])

