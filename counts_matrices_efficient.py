import anndata as ad
import pandas as pd
import numpy as np
import os
import subprocess
import time

# DEPRECATED

# Parameters
size = 1000

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

# Extract boolean masks for RNA and ATAC features and 0 count columns
rna_mask = adata.var["feature_types"] == "GEX"
atac_mask = adata.var["feature_types"] == "ATAC"
zero_counts = adata.layers["counts"].sum(axis=0) != 0
print("Shape of zero_counts filtering:")
print(adata[:, zero_counts].layers["counts"].shape)
zero_counts = np.array(zero_counts).flatten()


# Print shapes of masks
print("RNA mask shape:")
print(rna_mask.shape)
print("ATAC mask shape:")
print(atac_mask.shape)
print("Zero counts shape:")
print(zero_counts.shape)

# Combine feature mask with zero count mask
rna_mask = rna_mask & zero_counts
atac_mask = atac_mask & zero_counts

# Extract RNA UMI count matrices
print('Extracting RNA UMI count matrices')
rna_umi_counts = adata[:, rna_mask].layers["counts"][:size, :].T
print("Shape of rna_umi_counts:")
print(rna_umi_counts.shape)

# Extract RNA index
rna_index = adata[:, rna_mask].var.index

# Convert sparse matrix to dense
print('Starting sparse matrix conversion at', time.ctime())
rna_umi_counts = rna_umi_counts.toarray()
print('Finished sparse matrix conversion at', time.ctime())

# Write to file
print('Starting RNA write at', time.ctime())
rna_df = pd.DataFrame(rna_umi_counts, index=rna_index)
rna_df.to_csv(os.path.join(OUTPUT_DIR, "RNA_countmatrix_" + str(size) + "cell_filtered.txt"), sep='\t', header=False)
print('Finished RNA write at', time.ctime())

# Extract ATAC Peak count matrices
print('Extracting ATAC Peak count matrices')
atac_peak_counts = adata[:, atac_mask].layers["counts"][:size, :].T
print("Shape of atac_peak_counts:")
print(atac_peak_counts.shape)

# Extract ATAC index
atac_index = adata[:, atac_mask].var.index

# Convert sparse matrix to dense
print('Starting sparse matrix conversion at', time.ctime())
atac_peak_counts = atac_peak_counts.toarray()
print('Finished sparse matrix conversion at', time.ctime())


print('Starting ATAC write at', time.ctime())
atac_df = pd.DataFrame(atac_peak_counts, index=atac_index)
atac_df.to_csv(os.path.join(OUTPUT_DIR, "ATAC_countmatrix_" + str(size) + "cell_filtered.txt"), sep='\t', header=False)
print('Finished ATAC write at', time.ctime())

# Extract cell types
cell_types = pd.read_csv('/home/mzr19001/outputs/cell_types.txt', header=None)

#Output cell types as a single column text file without identifiers
cell_types[:size].to_csv(os.path.join(OUTPUT_DIR, "cell_types_" + str(size) + ".txt"), sep='\t', header=False, index=False)

#Ending
print('Finished after ', time.time() - start)