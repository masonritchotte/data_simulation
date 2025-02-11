import pandas as pd

def calculate_sparsity(file_path):
    df = pd.read_csv(file_path, delimiter='\t')  # Adjust delimiter as needed
    total_elements = df.size
    zero_elements = (df == 0).sum().sum()
    sparsity_percentage = (zero_elements / total_elements) * 100
    return sparsity_percentage

file_path = '/home/mzr19001/outputs/RNA_countmatrix_2500cell_filtered.txt'  # Replace with your file path
sparsity = calculate_sparsity(file_path)
print(f"Sparsity Percentage: {sparsity:.2f}%")