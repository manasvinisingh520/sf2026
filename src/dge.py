"""
Differential Gene Expression (DGE) Analysis Class

This module provides a class-based interface for performing DGE analysis
using PyDESeq2 (DESeq2 implementation in Python).
"""

import pandas as pd
import numpy as np
from typing import Optional, Tuple, Dict, List
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
    
    Attributes:
    -----------
    method : str
        DGE method to use (default: "deseq2")
    data_file : AnnData
        AnnData object containing gene expression data
    threshold : dict
        Thresholds for significance filtering. Default: 
        {'padj': 0.05, 'log2fc': 1.0}
    topK : int
        Number of top genes to return (default: 10)
    output_file : str
        Output file name for results (default: "dge_results.csv")
    condition : str
        Column name in adata.obs for the condition/group variable
    group1 : str
        Name of the first group to compare (numerator in fold change)
    group2 : str
        Name of the second group to compare (denominator in fold change)
    min_counts : int
        Minimum total counts across samples to keep a gene (default: 10)
    n_cpus : int
        Number of CPUs to use for parallel processing (default: 4)
    control_genes : Optional[List[str]]
        List of control/housekeeping genes for normalization (default: None)
    """
    
    def __init__(
        self,
        data_file: AnnData,
        condition: str,
        group1: str,
        group2: str,
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
        group1 : str
            Name of group 1 (numerator in fold change calculation)
            Example: "treated", "disease", "B"
        group2 : str
            Name of group 2 (denominator in fold change calculation)
            Example: "control", "healthy", "A"
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
        
        # Store parameters
        self.data_file = data_file
        self.condition = condition
        self.group1 = group1
        self.group2 = group2
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
        
        # Validate condition column exists
        if condition not in data_file.obs.columns:
            raise ValueError(
                f"Condition column '{condition}' not found in adata.obs. "
                f"Available columns: {list(data_file.obs.columns)}"
            )
        
        # Validate groups exist
        unique_conditions = data_file.obs[condition].unique()
        if group1 not in unique_conditions:
            raise ValueError(
                f"Group '{group1}' not found in condition '{condition}'. "
                f"Available groups: {list(unique_conditions)}"
            )
        if group2 not in unique_conditions:
            raise ValueError(
                f"Group '{group2}' not found in condition '{condition}'. "
                f"Available groups: {list(unique_conditions)}"
            )
        
        # Initialize results attributes
        self.dds = None  # DeseqDataSet object
        self.stat_res = None  # DeseqStats object
        self.results = None  # Results DataFrame
        self.significant_genes = None  # Filtered significant genes
        
        # Store extracted counts and metadata
        self.counts_df = None
        self.metadata_df = None
    
    def _prepare_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Prepare count matrix and metadata from AnnData object.
        
        Returns:
        --------
        counts_df : pd.DataFrame
            Count matrix with genes as rows, samples as columns
        metadata_df : pd.DataFrame
            Sample metadata with condition information
        """
        print("Preparing data from AnnData object...")
        
        # Extract counts matrix (cells × genes)
        # Convert to genes × cells (samples) format for DESeq2
        if hasattr(self.data_file.X, 'toarray'):
            # Sparse matrix
            counts_matrix = self.data_file.X.toarray().T  # Transpose: genes × cells
        else:
            # Dense matrix
            counts_matrix = self.data_file.X.T  # Transpose: genes × cells
        
        # Create counts DataFrame with gene names as index
        self.counts_df = pd.DataFrame(
            counts_matrix,
            index=self.data_file.var_names,
            columns=self.data_file.obs_names
        )
        
        # Extract metadata
        self.metadata_df = self.data_file.obs.copy()
        
        # Filter lowly expressed genes
        if self.min_counts > 0:
            gene_sums = self.counts_df.sum(axis=1)
            n_genes_before = len(self.counts_df)
            self.counts_df = self.counts_df[gene_sums >= self.min_counts]
            n_genes_after = len(self.counts_df)
            print(f"Filtered genes: {n_genes_before} -> {n_genes_after} "
                  f"(kept genes with >= {self.min_counts} total counts)")
        
        print(f"Counts matrix shape: {self.counts_df.shape} (genes × cells)")
        print(f"Metadata shape: {self.metadata_df.shape}")
        print(f"Condition '{self.condition}': {self.metadata_df[self.condition].value_counts().to_dict()}")
        
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
        
        # Extract results
        print(f"\nExtracting results: {self.group1} vs {self.group2}")
        self.stat_res = DeseqStats(
            self.dds,
            contrast=[self.condition, self.group1, self.group2],
            inference=inference,
            n_cpus=self.n_cpus
        )
        
        # Generate summary
        self.stat_res.summary()
        
        # Get results DataFrame
        self.results = self.stat_res.results_df
        
        print(f"\nDGE analysis complete. Total genes tested: {len(self.results)}")
        
        return self.results
    
    def find_significant_genes(self) -> pd.DataFrame:
        """
        Find significant differentially expressed genes based on thresholds.
        
        Returns:
        --------
        significant_genes : pd.DataFrame
            DataFrame with significant genes filtered by padj and log2fc thresholds
        """
        if self.results is None:
            raise ValueError("Results not available. Run run_dge() first.")
        
        padj_threshold = self.threshold.get('padj', 0.05)
        log2fc_threshold = self.threshold.get('log2fc', 1.0)
        
        print(f"\nFiltering significant genes...")
        print(f"  Thresholds: padj < {padj_threshold}, |log2FC| > {log2fc_threshold}")
        
        # Filter significant genes
        mask = (
            (self.results['padj'] < padj_threshold) &  # Statistically significant
            (abs(self.results['log2FoldChange']) > log2fc_threshold)  # Biologically significant
        )
        
        self.significant_genes = self.results[mask].copy()
        
        print(f"  Significant DEGs found: {len(self.significant_genes)}")
        
        # Print summary by direction
        if len(self.significant_genes) > 0:
            upregulated = (self.significant_genes['log2FoldChange'] > 0).sum()
            downregulated = (self.significant_genes['log2FoldChange'] < 0).sum()
            print(f"    Upregulated ({self.group1} > {self.group2}): {upregulated}")
            print(f"    Downregulated ({self.group1} < {self.group2}): {downregulated}")
        
        return self.significant_genes
    
    def get_top_genes(self, k: Optional[int] = None, sort_by: str = 'padj') -> pd.DataFrame:
        """
        Get top K genes sorted by significance.
        
        Parameters:
        -----------
        k : int, optional
            Number of top genes to return. If None, uses self.topK (default: None)
        sort_by : str, optional
            Column to sort by. Options: 'padj', 'pvalue', 'log2FoldChange' (default: 'padj')
        
        Returns:
        --------
        top_genes : pd.DataFrame
            Top K genes sorted by specified column
        """
        if self.significant_genes is None:
            # Use all results if significant_genes not filtered yet
            if self.results is None:
                raise ValueError("Results not available. Run run_dge() first.")
            genes_df = self.results
        else:
            genes_df = self.significant_genes
        
        if k is None:
            k = self.topK
        
        # Validate sort_by column
        if sort_by not in genes_df.columns:
            raise ValueError(
                f"Sort column '{sort_by}' not found. "
                f"Available columns: {list(genes_df.columns)}"
            )
        
        # Sort and get top K
        ascending = True if sort_by in ['padj', 'pvalue'] else False
        
        top_genes = genes_df.sort_values(sort_by, ascending=ascending).head(k)
        
        print(f"\nTop {k} genes sorted by {sort_by}:")
        print(top_genes[['log2FoldChange', 'padj', 'baseMean']].head(k))
        
        return top_genes
    
    def save_results(self, filename: Optional[str] = None) -> str:
        """
        Save DGE results to CSV file.
        
        Parameters:
        -----------
        filename : str, optional
            Output filename. If None, uses self.output_file (default: None)
        
        Returns:
        --------
        filename : str
            Path to saved file
        """
        if self.results is None:
            raise ValueError("Results not available. Run run_dge() first.")
        
        if filename is None:
            filename = self.output_file
        
        self.results.to_csv(filename)
        print(f"\nResults saved to '{filename}'")
        
        return filename
    
    def run_full_analysis(self) -> Dict[str, pd.DataFrame]:
        """
        Run complete DGE analysis pipeline.
        
        This method runs all steps:
        1. Run DGE analysis
        2. Find significant genes
        3. Get top K genes
        4. Save results
        
        Returns:
        --------
        results_dict : dict
            Dictionary containing:
            - 'all_results': All DGE results
            - 'significant_genes': Significant DEGs
            - 'top_genes': Top K genes
        """
        print("=" * 60)
        print("Running Full DGE Analysis Pipeline")
        print("=" * 60)
        
        # Step 1: Run DGE
        all_results = self.run_dge()
        
        # Step 2: Find significant genes
        significant = self.find_significant_genes()
        
        # Step 3: Get top genes
        top_genes = self.get_top_genes()
        
        # Step 4: Save results
        self.save_results()
        
        results_dict = {
            'all_results': all_results,
            'significant_genes': significant,
            'top_genes': top_genes
        }
        
        print("\n" + "=" * 60)
        print("Analysis Complete!")
        print("=" * 60)
        
        return results_dict
    
    def __repr__(self):
        """String representation of the DGEAnalyzer object."""
        status = "Analysis run" if self.results is not None else "Not run"
        return (
            f"DGEAnalyzer(\n"
            f"  method='{self.method}',\n"
            f"  condition='{self.condition}',\n"
            f"  comparison='{self.group1} vs {self.group2}',\n"
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
    
    This example demonstrates how to use the class with an AnnData object.
    """
    # Example: Load AnnData object (replace with your actual data)
    # import scanpy as sc
    # adata = sc.read_h5ad("your_data.h5ad")
    
    # Or create a synthetic example
    from anndata import AnnData
    import numpy as np
    
    # Create synthetic data for demonstration
    n_cells = 100
    n_genes = 1000
    
    # Synthetic count matrix
    counts = np.random.negative_binomial(n=10, p=0.3, size=(n_cells, n_genes))
    
    # Create metadata with conditions
    metadata = pd.DataFrame({
        'condition': ['A'] * (n_cells // 2) + ['B'] * (n_cells // 2),
        'sample_id': [f'S{i}' for i in range(n_cells)]
    })
    
    # Create AnnData object
    adata = AnnData(counts)
    adata.obs = metadata
    adata.var_names = [f'Gene_{i}' for i in range(n_genes)]
    adata.obs_names = [f'Cell_{i}' for i in range(n_cells)]
    
    # Create DGE analyzer
    dge = DGEAnalyzer(
        data_file=adata,
        condition="condition",
        group1="B",
        group2="A",
        method="deseq2",
        threshold={'padj': 0.05, 'log2fc': 0.5},
        topK=10,
        output_file="dge_results.csv",
        min_counts=10,
        n_cpus=2
    )
    
    # Run full analysis
    results = dge.run_full_analysis()
    
    # Access results
    print("\nAccessing results:")
    print(f"Total genes tested: {len(results['all_results'])}")
    print(f"Significant DEGs: {len(results['significant_genes'])}")
    print(f"Top {dge.topK} genes:")
    print(results['top_genes'])

