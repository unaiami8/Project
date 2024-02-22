library(Matrix)
library(ggplot2)
library(limma)
library(anndata)
library(Seurat)

covid_adata = read_h5ad("/active/debruinz_project/gautam_subedi/adata_covid.h5ad")
normal_adata = read_h5ad("/active/debruinz_project/gautam_subedi/adata_normal.h5ad")
covid_exp_mat = t(covid_adata$X)
normal_exp_mat = t(normal_adata$X)
covid_sub <- covid_exp_mat[,sample(1:ncol(covid_exp_mat),300000)]
normmal_sub <- normal_exp_mat[,sample(1:ncol(normal_exp_mat),300000)]
combined_mat <- cbind(covid_sub, normal_sub)
A <- as(combined_mat, "CsparseMatrix")


#Create design matrix, fit limma and then result
# create design matrix, create voom, create contrast matrix, fit limma, result 
design <- model.matrix(~ 0 + factor(c(rep("COVID-19", ncol(combined_mat)), rep("normal", ncol(combined_mat))))
