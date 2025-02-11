import anndata as ad
import pandas as pd
import numpy as np
import os
import time

# Parameters
size = 2500  # number of cells to subset

# Input H5ad file and Output directory
H5AD_FILE = "/home/mzr19001/datasets/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
OUTPUT_DIR = "/home/mzr19001/outputs/"

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

print('Starting file read at', time.ctime())
start = time.time()

# Load the AnnData object
adata = ad.read_h5ad(H5AD_FILE)
print('Finished file read at', time.ctime())

# Define feature masks based on modality
rna_mask = adata.var["feature_types"] == "GEX"
atac_mask = adata.var["feature_types"] == "ATAC"

#################################
# Process RNA modality
#################################
print('Extracting RNA UMI count matrices')

# Subset the first `size` cells and RNA features
rna_matrix = adata[:size, rna_mask].layers["counts"]

# Convert the sparse matrix to a dense numpy array
rna_matrix_dense = rna_matrix.toarray()

# Compute variance for each gene (feature) across cells (axis=0)
rna_variance = rna_matrix_dense.var(axis=0)

# Select the indices for the top 1000 genes with highest variance
top_1000_idx = np.argsort(rna_variance)[-1000:][::-1]

# Subset the RNA count matrix to the top genes
rna_matrix_filtered = rna_matrix[:, top_1000_idx]

# Get the corresponding gene names
rna_index = adata[:, rna_mask].var.index[top_1000_idx]

# Transpose so that features (genes) become rows (shape: [# genes, size])
rna_umi_counts = rna_matrix_filtered.T.toarray()
print("Shape of filtered RNA UMI count matrix:", rna_umi_counts.shape)

# Write RNA data to file
print('Starting RNA write at', time.ctime())
rna_df = pd.DataFrame(rna_umi_counts, index=rna_index)
rna_filename = os.path.join(OUTPUT_DIR, f"RNA_countmatrix_{size}cell_top1000genes.txt")
rna_df.to_csv(rna_filename, sep='\t', header=False)
print('Finished RNA write at', time.ctime())

#################################
# Process ATAC modality
#################################
print('Extracting ATAC Peak count matrices')

# Subset the first `size` cells and ATAC features
atac_matrix = adata[:size, atac_mask].layers["counts"]

# Convert the sparse matrix to a dense numpy array
atac_matrix_dense = atac_matrix.toarray()

# Compute variance for each ATAC feature (peak) across cells (axis=0)
atac_variance = atac_matrix_dense.var(axis=0)

# Select the indices for the top 5000 peaks with highest variance
top_5000_idx = np.argsort(atac_variance)[-5000:][::-1]

# Subset the ATAC count matrix to the top peaks
atac_matrix_filtered = atac_matrix[:, top_5000_idx]

# Get the corresponding peak names
atac_index = adata[:, atac_mask].var.index[top_5000_idx]

# Transpose so that features (peaks) become rows (shape: [# peaks, size])
atac_peak_counts = atac_matrix_filtered.T.toarray()
print("Shape of filtered ATAC peak count matrix:", atac_peak_counts.shape)

# Write ATAC data to file
print('Starting ATAC write at', time.ctime())
atac_df = pd.DataFrame(atac_peak_counts, index=atac_index)
atac_filename = os.path.join(OUTPUT_DIR, f"ATAC_countmatrix_{size}cell_top5000peaks.txt")
atac_df.to_csv(atac_filename, sep='\t', header=False)
print('Finished ATAC write at', time.ctime())

#################################
# Process cell types
#################################
# Read in cell types from file and subset to the first `size` cells.
cell_types = pd.read_csv('/home/mzr19001/outputs/cell_types.txt', header=None)
cell_types_filename = os.path.join(OUTPUT_DIR, f"cell_types_{size}.txt")
cell_types[:size].to_csv(cell_types_filename, sep='\t', header=False, index=False)

# Finished
print('Finished after ', time.time() - start)