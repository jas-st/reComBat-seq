library(anndataR)
library(reComBatseq)

# read in the dataset
dat <- read_h5ad("tutorial/Adipose_tissue_exp.h5ad")
sce <- dat$as_SingleCellExperiment()
counts_df <- sce@assays@data@listData[["X"]]
colnames(counts_df) <- sce@colData@rownames

batches <- sce[["batch"]]
group <- sce[["disease"]]

# which batches contain more than 1 sample
valid_batches = names(table(batches)[table(batches) > 1])
keep_lst_batches = which(batches %in% valid_batches)
counts_df_reduced <- counts_df[,keep_lst_batches]

batches <- droplevels(batches[keep_lst_batches])
group <- droplevels(group[keep_lst_batches])

# batch correction
recombatseq_df <- reComBat_seq(counts_df_reduced, batch = batches, group = group,
                             lambda_reg=0.8, alpha_reg=0.3)


count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
batch <- c(rep(1, 4), rep(2, 4))

adjusted <- reComBat_seq(count_matrix, batch=batch, group=NULL)
