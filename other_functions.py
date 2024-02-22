# another approach, This can be revised to use for sparse matrix
gene_names = covid_adata.var.feature_name
num_genes = covid_exp_mat.shape[1]
gene_indices = np.arange(0, num_genes)
t_statistics = np.zeros(num_genes)
p_values = np.zeros(num_genes)
significant_gene_names = []
significant_t_statistics = []
significant_p_values = []
for gene_index in range(num_genes):
    covid_gene_data = covid_exp_mat[:, gene_index].data
    normal_gene_data = normal_exp_mat[:, gene_index].data
    if len(covid_gene_data) >= 2 and len(normal_gene_data) >= 2:
        t_statistic, p_value = ttest_ind(covid_gene_data, normal_gene_data, equal_var=False)
        t_statistics[gene_index] = t_statistic
        p_values[gene_index] = p_value
        if p_value < 0.05:
            significant_gene_names.append(gene_names.values[gene_index])
            significant_t_statistics.append(t_statistic)
            significant_p_values.append(p_value)
    else:
        t_statistics[gene_index] = np.nan
        p_values[gene_index] = np.nan