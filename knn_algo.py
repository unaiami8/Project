import scanpy as sc
import pandas as pd
import os
from sklearn.metrics import pairwise_distances

# Load the h5ad file
adata = sc.read('your_file.h5ad')

# Perform nearest neighbor analysis on the h5ad data
sc.pp.neighbors(adata)

# Directory containing CSV files
csv_dir = '/active/debruinz_project/CellCensus/Python'

# Initialize variables to store the most similar CSV file and its similarity score
most_similar_csv = None
best_similarity_score = float('inf')  # Set to infinity for minimizing distance

# Loop through CSV files in the directory
for file in os.listdir(csv_dir):
    if file.endswith('.csv'):
        # Load CSV file
        csv_data = pd.read_csv(os.path.join(csv_dir, file))
        
        # Convert to AnnData object (assuming CSV contains similar data structure)
        adata_csv = sc.AnnData(csv_data)
        
        # Perform nearest neighbor analysis on the CSV data
        sc.pp.neighbors(adata_csv)
        
        # Calculate similarity score between the neighborhood graphs
        # Here, we are using Euclidean distance between adjacency matrices
        similarity_score = pairwise_distances(adata.obsp['connectivities'], adata_csv.obsp['connectivities']).mean()
        
        # Update most similar CSV file if similarity score is better
        if similarity_score < best_similarity_score:
            most_similar_csv = file
            best_similarity_score = similarity_score

# Print the most similar CSV file and its similarity score
print("Most similar CSV file:", most_similar_csv)
print("Similarity score:", best_similarity_score)
