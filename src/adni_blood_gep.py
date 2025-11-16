import pandas as pd
import os

# Get the path to the CSV file
data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
csv_file = os.path.join(data_dir, 'ADNI_Gene_Expression_Profile.csv')

# Read first 9 rows without headers to examine the structure
df_temp = pd.read_csv(csv_file, header=None, nrows=9)

# Extract Subject IDs from row 3 (index 2) - these are in columns 3 onwards
subject_ids = df_temp.iloc[2, 3:].tolist()

# Create column names: first 3 are ProbeSet, LocusLink, Gene, then Subject IDs
column_names = ['ProbeSet', 'LocusLink', 'Gene'] + subject_ids

# Read the CSV file starting from row 10 (skip first 9 rows which includes the old header)
df = pd.read_csv(csv_file, skiprows=9, names=column_names)

# Print the new column names
print("Column Names:")
print(df.columns.tolist()[:10])  # Show first 10 column names
print(f"\nTotal columns: {len(df.columns)}")
print(f"\nDataFrame shape: {df.shape}")
print("\nFirst few rows:")
print(df.head())

# List and unique values of the Gene column
print("\nUnique values in Gene column:")
print(df['Gene'].unique())
# Also print length of the Gene column
print(f"\nLength of Gene column: {len(df['Gene'])}")

# Create list of known genes
my_genes_text = """APOE Smooth muscle cells, Meningeal fibroblasts
PICALM * Brain endothelium (arterial, capillary, and venous)
CLU Meningeal fibroblasts, Ependymal cells
ABCA7 T cells
PTK2B T cells
PLCG2 * Brain endothelium (arterial)
HLA-DRB1 * Brain endothelium (arterial)
CD2AP * Brain endothelium (arterial, capillary, and venous)
SLC24A4 Ependymal cells
RIN3 T cells
ADAMTS1 * Smooth muscle cells, Pericytes, Brain endothelium (arterial)
ADAMTS4 Smooth muscle cells, Pericytes
FERMT2 Smooth muscle cells, Pericytes
SCIMP Ependymal cells
CLNK T cell
ECHDC3 Perivascular fibroblasts
TNIP1 * Brain endothelium (capillary and venous), T cells
ABCA1 * Brain endothelium (venous), Perivascular fibroblasts
USP6NL * Brain endothelium (capillary and venous)
INPP5D * Brain endothelium (capillary), T cells
ACE * Brain endothelium (arterial and capillary)
IQCK Ependymal cells
ABI3 T cells
HESX1 Meningeal fibroblasts
FHL2 Perivascular fibroblasts
CHRNE Perivascular fibroblasts, T cells
AGRN Pericytes"""
# Parse the text into a list of genes
my_genes = my_genes_text.split('\n')
my_genes = [gene.split(' ')[0].strip() for gene in my_genes]
# Print the list of genes
print(my_genes)
print(f"\nLength of my_genes: {len(my_genes)}")

## Find 
my_genes_rows = df[df['Gene'].isin(my_genes)]
# Print the rows where the Gene column matches the list of known genes
print(my_genes_rows)
print(f"\nLength of my_genes_rows: {len(my_genes_rows)}")
