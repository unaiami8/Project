import cellxgene_census
import pandas as pd
import scanpy as sc
import math
import tiledb
import anndata
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import anndata2ri
import rpy2.robjects as ro

#Import R libraries
base = importr('base')
stats = importr('stats')
limma = importr('limma')
writexl = importr('writexl')

# interactive node
#srun -n 4 --time=100:00:00 --pty bash -i
#BIGMEM: srun -p bigmem -n 4 --time=100:00:00 --pty bash -i

# my directory
# cd active/debruinz_project/gautam_subedi

#opening soma
census = cellxgene_census.open_soma()

# organ = 'string'

# COVID-19 adata with filter
adata_covid = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    obs_value_filter= "disease == 'COVID-19'and tissue_general == 'lung' and is_primary_data == True" ,
    column_names = {"obs": ["assay", "cell_type", "tissue", "tissue_general", "disease"]},
)

#Normal adata with same filter
adata_normal = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    obs_value_filter= "disease == 'normal'and tissue_general == 'lung' and is_primary_data == True" ,
    column_names = {"obs": ["assay", "cell_type", "tissue", "tissue_general", "disease"]},
)

#Saving adata as h5ad file
output_path1 = "/active/debruinz_project/gautam_subedi/adata_covid.h5ad"
output_path2 = "/active/debruinz_project/gautam_subedi/adata_normal.h5ad"
adata_covid.write(output_path1)
adata_normal.write(output_path2)

#Read anndata
covid_adata = anndata.read_h5ad("/active/debruinz_project/gautam_subedi/adata_covid.h5ad")
normal_adata = anndata.read_h5ad("/active/debruinz_project/gautam_subedi/adata_normal.h5ad")

#Fetching expression matrix
covid_exp_mat = covid_adata.X
normal_exp_mat = normal_adata.X

#checking gene_names and cell_type of covid_data
covid_adata.var.feature_name
covid_adata.obs.cell_type

#checking gene_names and cell_type of normal_data and covid and non-covid have same feature length and ids
normal_adata.var.feature_name
normal_adata.obs.cell_type

#unique cell_type in covid_adata
adata_covid.obs['cell_type'].unique  #70 categories
adata_normal.obs.feature_name.unique  #60664

## unique cell_type in normal
adata_normal.obs['cell_type'].unique() #184 catrgoreis

# This section shows common expression matrix containing genes and cells that are present in both.

# CELL METADATA: Common cell_type
unique_cell_types_covid = set(adata_covid.obs['cell_type'].unique())
unique_cell_types_normal = set(adata_normal.obs['cell_type'].unique())
common_cell_types = unique_cell_types_covid.intersection(unique_cell_types_normal)
#print(common_cell_types)  = 66 common cell types

# GENE METADATA: Output shows that they have same gene name and also at same position with same dimension
covid_feature_names = set(covid_adata.var['feature_name'])
normal_feature_names = set(normal_adata.var['feature_name'])
common_feature_names = covid_feature_names.intersection(normal_feature_names)
#len(common_feature_names)
#60664

# Define the number of cells or genes you want to subsample
num_cells = 20000

# Randomly select a subset of cells
selected_cells_indices = np.random.choice(covid_exp_mat.shape[0], size=num_cells, replace=False)
subsampled_covid_exp_mat = covid_exp_mat[ selected_cells_indices, :]
selected_cells_indices = np.random.choice(normal_exp_mat.shape[0], size=num_cells, replace=False)
subsampled_noraml_exp_mat = normal_exp_mat[selected_cells_indices, :]


# Load limma library in R
ro.r("library(limma)")
ro.r("library(Matrix)")

# These are conversions from sparse matrices to R
mmwrite("normal_mat_R.mtx", subsampled_noraml_exp_mat)
mmwrite("covid_mat_R.mtx", subsampled_covid_exp_mat)

# Turning sparse matrices into the following variables in R
ro.r("covid_data <- readMM('covid_mat_R.mtx')")
ro.r("non_covid_data <- readMM('normal_mat_R.mtx')")

# Combine data into a data frame or matrix with columns corresponding to each group
# Assuming non_covid_data and covid_data are numeric vectors or matrices of equal length

# Assuming 'combined_data' is your variable or matrix
ro.r("covid_data <- t(covid_data)")
ro.r("non_covid_data <- t(non_covid_data)")


ro.r("combined_data <- cbind(non_covid_data, covid_data)")
ro.r("combined_data <- t(combined_data)")


ro.r("design <- cbind(NonCovid = rep(1, ncol(non_covid_data)), Covid = rep(0, ncol(covid_data)))")

ro.r("fit <- lmFit(combined_data, design)")

ro.r("fit <- eBayes(fit)")

ro.r("results <- decideTests(fit)")
ro.r("significant_genes <- rownames(combined_data)[results$FDR <= 0.05]")