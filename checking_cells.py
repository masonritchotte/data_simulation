import pandas as pd
import numpy as np

# Load the data
df = pd.read_csv('/home/mzr19001/outputs/ATAC_countmatrix_2500cell_filtered.txt', sep='\t', header=None, index_col=0)


# Check for zero count columns
zero_counts = df.sum(axis=1) == 0

print("Number of zero count columns:")
print(zero_counts.sum())
