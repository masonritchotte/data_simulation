import pandas as pd
import matplotlib.pyplot as plt

FILE_TAG = "2500.13.True"

# Specify the path to the cell types text file
CELL_TYPES_FILE = "/home/mzr19001/data_simulation/Results/synthetic/scMultiOmics.scDesign3Simulated.CellTypeLabel.{}.txt".format(FILE_TAG)
# Specify the path to save the plotted image
PLOT_SAVE_PATH = "/home/mzr19001/data_simulation/Distribution Images/cell_types_distribution.{}.png".format(FILE_TAG)

def main():
    # Read the file (assumes no header, tab-separated, single column)
    cell_types = pd.read_csv(CELL_TYPES_FILE, header=None, sep='\t').iloc[:,0]
    
    # Compute distribution
    distribution = cell_types.value_counts().sort_index()
    print("Cell types distribution:")
    print(distribution)
    
    # Plot histogram as a bar chart
    plt.figure(figsize=(10,6))
    plt.bar(distribution.index.astype(str), distribution.values)
    plt.xlabel("Cell Type")
    plt.ylabel("Count")
    plt.title("Distribution of Cell Types")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(PLOT_SAVE_PATH)  # Save the plot to the specified path
    print(f"Plot saved to {PLOT_SAVE_PATH}")
    
    # Text-based histogram
    # print("\nText-based Histogram:")
    # max_val = distribution.max()
    # max_bar_width = 50  # maximum width in characters
    
    # for cell_type, count in distribution.items():
    #     bar_length = int((count / max_val) * max_bar_width)
    #     print(f"{str(cell_type):>10} | {'*' * bar_length} ({count})")

if __name__ == "__main__":
    main()
