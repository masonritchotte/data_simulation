import anndata as ad
import pandas as pd
import numpy as np
import os
import time

# Parameters for cell sampling and feature selection
size = 5000                     # Total number of cells to sample
desired_num_cell_types = 4      # Number of unique cell types to include (complexity)
balanced = True                 # True: sample equal numbers from each selected cell type; False: imbalanced sampling

# Input H5ad file and Output directory
H5AD_FILE = "/home/mzr19001/datasets/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
OUTPUT_DIR = "/home/mzr19001/data_simulation/Results"
CELL_TYPES_FILE = "/home/mzr19001/data_simulation/cell_types.txt"  # Assumed one-column file without header

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

print('Starting file read at', time.ctime())
start = time.time()

# Load the AnnData object
adata = ad.read_h5ad(H5AD_FILE)
print('Finished file read at', time.ctime())

# Read cell types (assumes one column, no header, and that the order matches adata.obs)
cell_types_df = pd.read_csv(CELL_TYPES_FILE, header=None)
# Convert to a Series for easier processing:
cell_types = cell_types_df.iloc[:, 0]

# Determine unique cell types available in the dataset
unique_cell_types = cell_types.unique()
print("Total unique cell types in data:", len(unique_cell_types))

# Choose a subset of cell types if desired_num_cell_types is less than the available ones
if desired_num_cell_types < len(unique_cell_types):
    chosen_cell_types = np.random.choice(unique_cell_types, size=desired_num_cell_types, replace=False)
else:
    chosen_cell_types = unique_cell_types
print("Selected cell types for sampling:", chosen_cell_types)

# Get indices of cells whose type is among the chosen types
allowed_cell_indices = cell_types.index[cell_types.isin(chosen_cell_types)].tolist()

# Now sample cell indices according to the "balanced" flag
if balanced:
    # For balanced sampling, aim for an equal number of cells per cell type
    per_type_count = size // desired_num_cell_types
    sampled_indices = []
    for ct in chosen_cell_types:
        ct_indices = cell_types.index[cell_types == ct].tolist()
        if len(ct_indices) < per_type_count:
            # Not enough cells of this type: use all available
            sampled_indices.extend(ct_indices)
        else:
            sampled_indices.extend(np.random.choice(ct_indices, per_type_count, replace=False))
    # If the total number of sampled cells is less than size, add extra cells randomly from the allowed pool.
    if len(sampled_indices) < size:
        remaining = list(set(allowed_cell_indices) - set(sampled_indices))
        additional_needed = size - len(sampled_indices)
        if len(remaining) >= additional_needed:
            additional = np.random.choice(remaining, additional_needed, replace=False)
            sampled_indices.extend(additional)
        else:
            sampled_indices.extend(remaining)
    # In case we sampled too many (can happen if size is not evenly divisible), trim the list.
    sampled_indices = sampled_indices[:size]
else:
    # For imbalanced sampling, sample 'size' cells at random from the allowed indices.
    if len(allowed_cell_indices) < size:
        print("Warning: not enough cells available in the selected cell types. Sampling all available cells.")
        sampled_indices = allowed_cell_indices
    else:
        sampled_indices = np.random.choice(allowed_cell_indices, size, replace=False)

# (Optional) Sort sampled indices to preserve original order
sampled_indices = np.sort(sampled_indices)
print("Number of cells sampled:", len(sampled_indices))

# Define feature masks based on modality (these are not affected by cell sampling)
rna_mask = adata.var["feature_types"] == "GEX"
atac_mask = adata.var["feature_types"] == "ATAC"

#################################
# Process RNA modality with variance filtering (top 1000 genes)
#################################
print('Extracting RNA UMI count matrices')

# Subset the AnnData object using the sampled indices for cells and the RNA features
rna_matrix = adata[sampled_indices, rna_mask].layers["counts"]

# Convert the sparse matrix to a dense array to compute variance
rna_matrix_dense = rna_matrix.toarray()
# Compute variance for each gene (feature) across the sampled cells (axis=0)
rna_variance = rna_matrix_dense.var(axis=0)

# Select the indices for the top 1000 genes with highest variance
top_1000_idx = np.argsort(rna_variance)[-1000:][::-1]

# Subset the RNA count matrix to these top genes
rna_matrix_filtered = rna_matrix[:, top_1000_idx]

# Get the corresponding gene names (from adata.var)
rna_index = adata[:, rna_mask].var.index[top_1000_idx]

# Transpose so that genes become rows (shape: [# genes, number of sampled cells])
rna_umi_counts = rna_matrix_filtered.T.toarray()
print("Shape of filtered RNA UMI count matrix:", rna_umi_counts.shape)

# Write RNA data to file
print('Starting RNA write at', time.ctime())
rna_df = pd.DataFrame(rna_umi_counts, index=rna_index)
rna_filename = os.path.join(OUTPUT_DIR, f"RNA_countmatrix_{size}cell_{desired_num_cell_types}_{balanced}_top1000genes.txt")
rna_df.to_csv(rna_filename, sep='\t', header=False)
print('Finished RNA write at', time.ctime())

#################################
# Process ATAC modality with variance filtering (top 5000 peaks)
#################################
print('Extracting ATAC Peak count matrices')

# Subset the AnnData object using the sampled cell indices for ATAC features
atac_matrix = adata[sampled_indices, atac_mask].layers["counts"]

# Convert the sparse matrix to a dense array for variance computation
atac_matrix_dense = atac_matrix.toarray()
# Compute variance for each ATAC peak (feature) across cells (axis=0)
atac_variance = atac_matrix_dense.var(axis=0)

# Select the indices for the top 5000 peaks with highest variance
top_5000_idx = np.argsort(atac_variance)[-5000:][::-1]

# Subset the ATAC count matrix to these top peaks
atac_matrix_filtered = atac_matrix[:, top_5000_idx]

# Get the corresponding peak names (from adata.var)
atac_index = adata[:, atac_mask].var.index[top_5000_idx]

# Transpose so that peaks become rows (shape: [# peaks, number of sampled cells])
atac_peak_counts = atac_matrix_filtered.T.toarray()
print("Shape of filtered ATAC peak count matrix:", atac_peak_counts.shape)

# Write ATAC data to file
print('Starting ATAC write at', time.ctime())
atac_df = pd.DataFrame(atac_peak_counts, index=atac_index)
atac_filename = os.path.join(OUTPUT_DIR, f"ATAC_countmatrix_{size}cell_{desired_num_cell_types}_{balanced}_top5000peaks.txt")
atac_df.to_csv(atac_filename, sep='\t', header=False)
print('Finished ATAC write at', time.ctime())

#################################
# Process cell types output
#################################
# Subset the cell types Series to the sampled indices (if you want to preserve the cell-type labels for the sampled cells)
sampled_cell_types = cell_types.iloc[sampled_indices]
cell_types_filename = os.path.join(OUTPUT_DIR, f"cell_types_{size}_{desired_num_cell_types}_{balanced}.txt")
sampled_cell_types.to_csv(cell_types_filename, sep='\t', header=False, index=False)

# Ending
print('Finished after ', time.time() - start)