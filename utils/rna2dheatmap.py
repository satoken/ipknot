import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def parse_rna_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = [line.strip() for line in lines if not line.startswith('#') and line.strip()]
    n = len(data)
    heatmap_matrix = np.zeros((n, n))

    for line in data:
        parts = line.split()
        base_index = int(parts[0]) - 1  # Convert to 0-based index
        for pair in parts[2:]:
            pair_index, probability = pair.split(':')
            pair_index = int(pair_index) - 1  # Convert to 0-based index
            probability = float(probability)
            heatmap_matrix[base_index, pair_index] = probability
            heatmap_matrix[pair_index, base_index] = probability  # Symmetric matrix

    return heatmap_matrix

def plot_heatmap(matrix, output_file):
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, cmap="viridis", cbar=True, square=True)
    plt.title("RNA Secondary Structure Heatmap")
    plt.xlabel("Base Index")
    plt.ylabel("Base Index")
    plt.savefig(output_file, format='png')
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rna_heatmap.py <input_file_path> <output_file_path>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    # データを解析してヒートマップをプロット
    heatmap_matrix = parse_rna_data(input_file_path)
    plot_heatmap(heatmap_matrix, output_file_path)

