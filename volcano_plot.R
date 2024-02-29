library(Matrix)
library(ggplot2)
library(limma)
library(anndata)
library(Seurat)

# Read the AnnData files
covid_adata <- read_h5ad("/active/debruinz_project/gautam_subedi/adata_covid.h5ad")
normal_adata <- read_h5ad("/active/debruinz_project/gautam_subedi/adata_normal.h5ad")

# Transpose Expression Matrices
covid_exp_mat <- t(covid_adata$X)
normal_exp_mat <- t(normal_adata$X)

# Optionally Subset Expression Matrices
covid_sub <- covid_exp_mat[, sample(1:ncol(covid_exp_mat), 20000)]
normal_sub <- normal_exp_mat[, sample(1:ncol(normal_exp_mat), 20000)]

# Combine Expression Matrices
combined_mat <- cbind(covid_sub, normal_sub)

# Convert to Sparse Matrix
A <- as(combined_mat, "CsparseMatrix")

# Create Design Matrix
design <- model.matrix(~ 0 + factor(c(rep("COVID-19", ncol(covid_sub)), rep("normal", ncol(normal_sub)))))

##########################################################################################################
##########################################################################################################
##########################################################################################################
# Start here:

# Load A (sparse matrix)
A <- readRDS("A_sparse_matrix.rds")

# Load design matrix
design <- readRDS("design_matrix.rds")

fit <- lmFit(A, design)

# Empirical Bayes Moderated t-Statistics
ebayes_fit <- eBayes(fit)

# Get p-values
results <- topTable(ebayes_fit, n = Inf)

# Change column names
colnames(results)[colnames(results) == "factor.c.rep..COVID.19...ncol.covid_sub....rep..normal...ncol.normal_sub....COVID.19"] <- "logFC"

# Volcano Plot
volcano_plot <- EnhancedVolcano(
  results,
  title = "Volcano Plot for Differential Expression Analysis",
  lab = rownames(results),   # Assuming gene names are in the row names
  x = "logFC",          # Specify the column name for log fold change
  y = "P.Value",        # Specify the column name for p-values
  xlim = c(-2, 2),      # Customize x-axis limits as needed
  ylim = c(0, 10),      # Customize y-axis limits as needed
  pCutoff = 0.05,       # P-value cutoff for highlighting
  FCcutoff = 1          # Log fold change cutoff for highlighting
)
