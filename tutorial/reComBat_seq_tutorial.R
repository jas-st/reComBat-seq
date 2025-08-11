library(anndataR)
library(reComBatseq)
library(SingleCellExperiment)
library(Rtsne)
library(ggplot2)

# read in the dataset
dat <- read_h5ad("tutorial/Esophagus_exp.h5ad")
sce <- dat$as_SingleCellExperiment()
counts_df <- assay(sce)
colnames(counts_df) <- sce@colData@rownames

## Remove genes with near 0 counts
keep_lst_genes <- which(apply(counts_df, 1, function(x){sum(x)>ncol(counts_df)}))
rm_genes <- setdiff(1:nrow(counts_df), keep_lst_genes)
counts_df_reduced <- counts_df[keep_lst_genes, ]

cat("Gene Amount after Removal: ", nrow(counts_df_reduced), "\n")

batches <- sce[["batch"]]
group <- sce[["disease"]]

# which batches contain more than 1 sample
valid_batches = names(table(batches)[table(batches) > 1])
keep_lst_batches = which(batches %in% valid_batches)
rm_batches <- setdiff(1:ncol(counts_df_reduced), keep_lst_batches)
counts_df_reduced <- counts_df_reduced[,keep_lst_batches]

batches <- droplevels(batches[keep_lst_batches])
group <- droplevels(group[keep_lst_batches])

# batch correction
recombatseq_df <- reComBat_seq(counts_df_reduced, batch = batches, group = group,
                             lambda_reg=0.8, alpha_reg=0.3)

# batch correction - confounded design
covmatdf <- as.data.frame(cbind(group, group))
covmatdf[] <- lapply(covmatdf, as.factor)
colnames(covmatdf) <- c("test1", "test2")
covmat_model <- model.matrix(~., data = covmatdf)[, -1]

recombatseq_df_conf <- reComBat_seq(counts_df_reduced, batch = batches, covar_mod = covmat_model,
                                   lambda_reg=0.8, alpha_reg=0.3)



# PLOT PCA - raw data
# Calculate tSNE using Rtsne(0 function)
tsne_raw <- Rtsne(t(counts_df_reduced), check_duplicates = FALSE,
                  perplexity = floor((ncol(counts_df_reduced) - 1) / 3))


# Conversion of matrix to dataframe
tsne_plot <- data.frame(x = tsne_raw$Y[,1],
                        y = tsne_raw$Y[,2],
                        batches = as.factor(batches),
                        group = as.factor(group))

# Plotting the plot using ggplot() function
ggplot2::ggplot(tsne_plot) + geom_point(aes(x=x, y=y, colour=group, shape=batches), size=3) +
  theme_bw()

ggsave("tutorial/PCA_raw.png", height=5, width=6, dpi=300)


# PLOT PCA - corrected data
# Calculate tSNE using Rtsne(0 function)
tsne_out <- Rtsne(t(recombatseq_df), check_duplicates = FALSE,
                  perplexity = floor((ncol(recombatseq_df) - 1) / 3))


# Conversion of matrix to dataframe
tsne_plot <- data.frame(x = tsne_out$Y[,1],
                        y = tsne_out$Y[,2],
                        batches = as.factor(batches),
                        group = as.factor(group))

# Plotting the plot using ggplot() function
ggplot2::ggplot(tsne_plot) + geom_point(aes(x=x, y=y, colour=group, shape=batches), size=3) +
  theme_bw()
ggsave("tutorial/PCA_corrected.png", height=5, width=6, dpi=300)


# PLOT PCA - corrected confounded data
# Calculate tSNE using Rtsne(0 function)
tsne_cor <- Rtsne(t(recombatseq_df_conf), check_duplicates = FALSE,
                  perplexity = floor((ncol(recombatseq_df_conf) - 1) / 3))


# Conversion of matrix to dataframe
tsne_plot <- data.frame(x = tsne_cor$Y[,1],
                        y = tsne_cor$Y[,2],
                        batches = as.factor(batches),
                        group = as.factor(group))

# Plotting the plot using ggplot() function
ggplot2::ggplot(tsne_plot) + geom_point(aes(x=x, y=y, colour=group, shape=batches), size=3) +
  theme_bw()

ggsave("tutorial/PCA_corrected_confounded.png", height=5, width=6, dpi=300)