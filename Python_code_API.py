import cellxgene_census
import pandas as pd
import scanpy as sc
import math
import tiledb
import anndata
from scipy.stats import ttest_ind
import numpy as np
import rpy2.robjects as ro
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import anndata2ri
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


#LIMMA, ANOVA
gene_names = covid_adata.var.feature_name
covid_gene_data = covid_exp_mat.toarray()
normal_gene_data = normal_exp_mat.toarray()
t_statistics, p_values = ttest_ind(covid_gene_data, normal_gene_data, equal_var=False)
significant_genes_value = p_values < 0.05
significant_gene_names = gene_names.values[significant_genes_value]
significant_t_statistics = t_statistics[significant_genes_value]
significant_p_values = p_values[significant_genes_value]
significant_gene_names =significant_gene_names.astype(str) #Not necessary if saving in dataframe
output_filename = "/active/debruinz_project/gautam_subedi/t_test_results.csv"
df = pd.DataFrame(data)
df.to_csv(output_filename, index=False)


# Convert pandas dataframes to R dataframes
pandas2ri.activate()

def limma_t_test(covid_exp_mat, normal_exp_mat):
    # Convert Python pandas dataframes to R dataframes
    covid_exp_mat_r = pandas2ri.py2ri(covid_exp_mat)
    normal_exp_mat_r = pandas2ri.py2ri(normal_exp_mat)

    # Run limma t-test
    design_matrix = robjects.r.cbind(robjects.IntVector([1]*covid_exp_mat_r.ncol), robjects.IntVector([0]*normal_exp_mat_r.ncol))
    fit = limma.lmFit(covid_exp_mat_r, design_matrix)
    contrast_matrix = robjects.r.matrix(robjects.FloatVector([1, -1]), nrow=1)
    contrast = robjects.r.contrasts.fit(fit, contrast_matrix)
    fit_eb = limma.eBayes(contrast)
    top_table = limma.topTable(fit_eb, coef=1, number=covid_exp_mat_r.nrow)
    
    # Extract results
    results = pandas2ri.ri2py_dataframe(top_table)
    
    return results

# Example usage:
# Assuming you have two pandas DataFrames: covid_exp_mat and normal_exp_mat
# results = limma_t_test(covid_exp_mat, normal_exp_mat)

## We can't convert the sparse matrices into R matrces for further analysis. WHen trying to make it 
## densa matrix and then to a exp_mat in R it doesnt load. 


## WHY DOESN'T THIS WORK ##
# Convert CSR matrices to dense arrays
covid_data_dense = covid_exp_mat.toarray()
non_covid_data_dense = normal_exp_mat.toarray()

# Transfer data to R
ro.numpy2ri.activate()  # Activate conversion between NumPy arrays and R arrays
r_covid_data = ro.r.matrix(covid_data_dense)
r_non_covid_data = ro.r.matrix(non_covid_data_dense)

# Load limma library in R
ro.r("library(limma)")

# Perform t-test using limma
ro.r("result <- eBayes(lmFit(cbind(covid_data, non_covid_data) ~ condition))")

sc.tl.rank_genes_groups(combined_adata, groups=['covid_exp_mat', 'normal_exp_mat'], method='ttest')  # Replace 'ttest' with desired test
