import pandas as pd
import numpy as np

FILE_TAG = "2500_13_True"

# Specify the path to the cell types text file
CELL_TYPES_FILE = "/home/mzr19001/data_simulation/Results/cell_types_{}.txt".format(FILE_TAG)

# Read the cell types file into a DataFrame
cell_types_df = pd.read_csv(CELL_TYPES_FILE, header=None, names=["CellType"])

# Count the occurrences of each cell type
cell_type_counts = cell_types_df["CellType"].value_counts()

# Sort the cell types
cell_type_counts = cell_type_counts.sort_index()

# Print the cell type counts
print(cell_type_counts)