import pandas as pd
import numpy as np
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# 1. Load Data (Using a synthetic example dataset for demonstration)
# In a real scenario, you would load your own 'counts.csv' and 'sample_info.csv'
from pydeseq2.utils import load_example_data
counts_df = load_example_data(modality="raw_counts", dataset="synthetic")
metadata = load_example_data(modality="metadata", dataset="synthetic")
print ("metadata: ", metadata)

# Optional: Filter out lowly expressed genes (e.g., keeping genes with at least 10 counts across all samples)
counts_df = counts_df[counts_df.sum(axis=1) >= 10]

print(f"Loaded counts data shape: {counts_df.shape}")
print(f"Loaded metadata shape: {metadata.shape}")

# 2. Create a DeseqDataSet object
# The design should match a column name in your metadata
# Use formulaic formula format: "~condition" (the ~ is optional but recommended)
design = "~condition"
inference = DefaultInference()

# OPTIONAL: Control genes for size factor calculation
# Control genes are used ONLY for size factor normalization (median-of-ratios method).
# They are genes expected to be invariant across conditions (e.g., housekeeping genes).
# How control genes are used:
#   1. When control_genes=None (default): ALL genes are used for size factor calculation
#   2. When control_genes is specified: ONLY those genes are used for size factor calculation
#   
# Why use control genes?
#   - Useful when many genes are differentially expressed, making median-of-ratios unreliable
#   - Use genes known to be stable across conditions (housekeeping genes like GAPDH, ACTB, etc.)
#   - Improves normalization robustness in cases with global expression changes
#
# Example: control_genes=['GAPDH', 'ACTB', 'RPL13A']  # If your counts_df has these gene names

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design=design,  # Use design instead of deprecated design_factors
    # ref_level is deprecated and no longer needed
    # Reference level is set in the contrast when calling DeseqStats
    # control_genes=None,  # Default: use all genes. Uncomment and specify list to use control genes
    inference=inference,
    n_cpus=4 # Adjust based on your system
)

# 3. Run the DGE analysis pipeline
print("Running DESeq2 pipeline...")
dds.deseq2()
print("Pipeline complete.")

# 4. Extract results
# This performs the Wald test for differential expression
stat_res = DeseqStats(
    dds,
    contrast=["condition", "B", "A"], # Comparing group 'B' vs 'A'
    inference=inference,
    n_cpus=4
)
stat_res.summary()
results = stat_res.results_df

# 5. Review and Filter Results
# 
# Understanding the key statistics in results:
#
# **log2FoldChange**: 
#   - The log2 (base 2) fold change in expression between conditions
#   - Formula: log2(mean_condition_B / mean_condition_A)
#   - Positive values: Gene is UP-regulated in condition B vs A (higher in B)
#   - Negative values: Gene is DOWN-regulated in condition B vs A (lower in B)
#   - Example: log2FoldChange = 2 means 2^2 = 4-fold increase in condition B
#   - Example: log2FoldChange = -1 means 2^(-1) = 0.5-fold (2x decrease) in condition B
#   - Threshold: |log2FoldChange| > 1 means |fold change| > 2x (doubling or halving)
#
# **padj** (adjusted p-value):
#   - Also called FDR (False Discovery Rate) or BH-adjusted p-value
#   - Statistical significance after correcting for multiple testing
#   - Controls the expected proportion of false positives among significant results
#   - Example: padj < 0.05 means < 5% false discovery rate
#   - More conservative than raw p-value (accounts for testing thousands of genes)
#   - Uses Benjamini-Hochberg procedure to adjust p-values
#
# **pvalue** (raw p-value):
#   - The probability of observing the data (or more extreme) if the null hypothesis is true
#   - Null hypothesis (H0): log2FoldChange = 0 (no differential expression)
#   - Alternative hypothesis (H1): log2FoldChange ≠ 0 (differential expression exists)
#   - Computed using the Wald test on the log2FoldChange estimate
#   - Formula: pvalue = P(|stat| > |observed_stat|) where stat ~ Normal(0, 1)
#   - Wald statistic: stat = log2FoldChange / lfcSE (follows standard normal distribution)
#   - pvalue < 0.05: Reject null hypothesis (gene is differentially expressed)
#   - IMPORTANT: Raw p-value does NOT account for multiple testing!
#     When testing thousands of genes, many will have pvalue < 0.05 by chance
#     This is why we use padj (adjusted p-value) instead
#
# How pvalue is computed in DESeq2/pydeseq2 (Wald test):
#   1. Estimate log2FoldChange (β) and its standard error (SE) from the GLM model
#   2. Compute Wald statistic: W = β / SE
#   3. Under null hypothesis (β = 0), W follows standard normal distribution N(0,1)
#   4. Calculate p-value: pvalue = 2 * P(Z > |W|) where Z ~ N(0,1)
#      (Two-tailed test: test if β ≠ 0)
#   5. The model accounts for:
#      - Negative binomial distribution of counts
#      - Dispersion estimates (gene-specific variance)
#      - Size factor normalization
#      - Design matrix (experimental conditions)
#
# **Other common columns in results**:
#   - baseMean: Average normalized count across all samples (used for filtering)
#   - lfcSE: Standard error of log2FoldChange estimate (uncertainty in FC)
#   - stat: Wald statistic (log2FoldChange / lfcSE)
#   - padj: Adjusted p-value (Benjamini-Hochberg FDR correction of pvalue)
#
# Filtering criteria:
#   - padj < 0.05: Statistically significant after multiple testing correction
#   - |log2FoldChange| > 1: Biologically significant (at least 2-fold change)

padj_threshold = 0.05  # FDR threshold (False Discovery Rate < 5%)
log2fc_threshold = 0.3   # Minimum |log2FoldChange| (i.e., |fold change| > 2x)

significant_genes = results[
    (results['padj'] < padj_threshold) &  # Statistically significant
    (abs(results['log2FoldChange']) > log2fc_threshold)  # Biologically significant
]

print(f"\nTotal genes tested: {results.shape[0]}")
print(f"Significant DEGs found (padj < {padj_threshold}, |log2FC| > {log2fc_threshold}): {significant_genes.shape[0]}")

# Display top significant genes
print("\nTop 10 significant genes:")
print(significant_genes.sort_values('padj').head(10))

# 6. Save results to a CSV file
results.to_csv("dge_results.csv")
print("\nResults saved to 'dge_results.csv'")
