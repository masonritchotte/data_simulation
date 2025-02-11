import anndata as ad
import pandas as pd
import numpy as np
import os
import subprocess
import time

# Parameters
size = 2500
# Input H5ad file and Output directory
H5AD_FILE = "/home/mzr19001/datasets/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
OUTPUT_DIR = "/home/mzr19001/outputs/"

# Ensure output directories exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

print('Starting file read at', time.ctime())
start = time.time()

# Load H5ad file
adata = ad.read_h5ad(H5AD_FILE)
print('Finished file read at', time.ctime())

# Extract boolean masks for RNA and ATAC features and 0 count columns (based on the full dataset)
rna_mask = adata.var["feature_types"] == "GEX"
atac_mask = adata.var["feature_types"] == "ATAC"
zero_counts = adata.layers["counts"].sum(axis=0) != 0
zero_counts = np.array(zero_counts).flatten()

# Combine feature mask with zero count mask (note: this is before subsetting cells)
rna_mask = rna_mask & zero_counts
atac_mask = atac_mask & zero_counts

#################################
# Process RNA modality
#################################
print('Extracting RNA UMI count matrices')

# Get the RNA count matrix for the first `size` cells
# (this matrix has shape: (size, # features))
rna_matrix = adata[:, rna_mask].layers["counts"][:size, :]

# Create a mask to filter out features with all zero counts in the selected cells.
# We sum along axis=0 (i.e. over cells) so each element corresponds to a feature.
nonzero_feature_mask = np.array(rna_matrix.sum(axis=0)).flatten() != 0

# Filter the RNA count matrix to only keep features with nonzero counts.
# Also update the RNA feature index accordingly.
rna_matrix_filtered = rna_matrix[:, nonzero_feature_mask]
rna_index = adata[:, rna_mask].var.index[nonzero_feature_mask]

# Transpose so that features become rows (shape: (# features, size))
# Then convert the filtered sparse matrix to a dense array.
rna_umi_counts = rna_matrix_filtered.T.toarray()

print("Shape of filtered rna_umi_counts:", rna_umi_counts.shape)

# Write RNA data to file
print('Starting RNA write at', time.ctime())
rna_df = pd.DataFrame(rna_umi_counts, index=rna_index)
rna_df.to_csv(os.path.join(OUTPUT_DIR, f"RNA_countmatrix_{size}cell_filtered.txt"), sep='\t', header=False)
print('Finished RNA write at', time.ctime())

#################################
# Process ATAC modality
#################################
print('Extracting ATAC Peak count matrices')

# Get the ATAC count matrix for the first `size` cells
atac_matrix = adata[:, atac_mask].layers["counts"][:size, :]

# Create a mask to filter out features with all zero counts in the selected cells.
nonzero_feature_mask = np.array(atac_matrix.sum(axis=0)).flatten() != 0

# Filter the ATAC count matrix and update the feature index accordingly.
atac_matrix_filtered = atac_matrix[:, nonzero_feature_mask]
atac_index = adata[:, atac_mask].var.index[nonzero_feature_mask]

# Transpose and convert to dense array.
atac_peak_counts = atac_matrix_filtered.T.toarray()

print("Shape of filtered atac_peak_counts:", atac_peak_counts.shape)

# Write ATAC data to file
print('Starting ATAC write at', time.ctime())
atac_df = pd.DataFrame(atac_peak_counts, index=atac_index)
atac_df.to_csv(os.path.join(OUTPUT_DIR, f"ATAC_countmatrix_{size}cell_filtered.txt"), sep='\t', header=False)
print('Finished ATAC write at', time.ctime())

#################################
# Process cell types
#################################
cell_types = pd.read_csv('/home/mzr19001/outputs/cell_types.txt', header=None)

# Output cell types as a single-column text file without identifiers (subset the first `size` cells)
cell_types[:size].to_csv(os.path.join(OUTPUT_DIR, f"cell_types_{size}.txt"), sep='\t', header=False, index=False)

# Ending
print('Finished after ', time.time() - start)
