#!/bin/bash

# Define job name
#SBATCH --job-name=limma_analysis

# Define output file
#SBATCH --output=limma_analysis.out

# Define error file
#SBATCH --error=limma_analysis.err

# Define walltime
#SBATCH --time=1:00:00

# Define partition
#SBATCH --partition=normal

# Load required modules
module load R

# Run R script
Rscript <<EOF
library(Matrix)
library(ggplot2)
library(limma)
library(anndata)
library(Seurat)

# Read the AnnData files
covid_adata <- 
read_h5ad("/active/debruinz_project/gautam_subedi/adata_covid.h5ad")
normal_adata <- 
read_h5ad("/active/debruinz_project/gautam_subedi/adata_normal.h5ad")

# Transpose Expression Matrices
covid_exp_mat <- t(covid_adata$X)
normal_exp_mat <- t(normal_adata$X)

# Optionally Subset Expression Matrices
covid_sub <- covid_exp_mat[, sample(1:ncol(covid_exp_mat), 300000)]
normal_sub <- normal_exp_mat[, sample(1:ncol(normal_exp_mat), 300000)]

# Combine Expression Matrices
combined_mat <- cbind(covid_sub, normal_sub)

# Convert to Sparse Matrix
A <- as(combined_mat, "CsparseMatrix")

# Create Design Matrix
design <- model.matrix(~ 0 + factor(c(rep("COVID-19", ncol(covid_sub)), 
rep("normal", ncol(normal_sub)))))
fit <- lmFit(A, design)
EOF

