# edgeR DE pipeline
edgeR_DEpipe <- function(cts, batch, group, include.batch, alpha.unadj, alpha.fdr, covar_incl=NULL, covar=NULL){
  y <- DGEList(counts=cts)
  y <- calcNormFactors(y, method="TMM")
  if(include.batch){
    cat("Including batch as covariate\n")
    if(is.null(group)){
      design <- model.matrix(~ as.factor(batch))
    }else{
      design <- model.matrix(~ as.factor(group) + as.factor(batch))
    }
  }else{
    cat("Default group as model matrix\n")
    if(is.null(group)){
      design <- model.matrix(~1, data=data.frame(t(cts)))
    }else{
      design <- model.matrix(~as.factor(group))
    }
  }
  if(!is.null(covar)){
    cat("Including surrogate variables or unwanted variation variables\n")
    design <- cbind(design, covar)
  }
  # include gender variable
  if(length(unique(covar_incl))>1){
    design <- cbind(design, Sex=model.matrix(~covar_incl)[,2])
  }

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef=2)
  de_res <- topTags(qlf, n=nrow(cts))$table

  de_called <- rownames(de_res)[de_res$PValue < alpha.unadj]
  de_called_fdr <- rownames(de_res)[de_res$FDR < alpha.fdr]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design))
}


##  DESeq2 pipeline
DESeq2_DEpipe <- function(counts_mat, batch, group, include.batch, alpha.unadj, alpha.fdr, covar_incl=NULL, covar=NULL){
  if(include.batch){
    cat("Including batch as covariate\n")
    col_data <- data.frame(Batch=as.factor(batch), Group=as.factor(group))
    design_formula <- ~Batch+Group
  }else if(!is.null(covar)){
    cat("Including surrogate variables or unwanted variation variables\n")
    col_data <- cbind(model.matrix(~as.factor(group)), covar)
    colnames(col_data)[2:ncol(col_data)] <- c("Group", paste0("Covar", 1:(ncol(col_data)-2)))
    rownames(col_data) <- colnames(counts_mat)
    col_data <- as.data.frame(col_data); col_data$Group <- as.factor(col_data$Group); col_data <- col_data[, -1]
    design_formula <- as.formula(paste("~", paste(colnames(col_data), collapse="+")))
  }else{
    cat("Default group as model matrix\n")
    col_data <- data.frame(Group=as.factor(group))
    design_formula <- ~Group
  }

  # deal with covar_incl
  if(!is.null(covar_incl)){
    covar_mod <- model.matrix(~factor(covar_incl))
    colnames(covar_mod)[-1] <- paste0("cov",1:(ncol(covar_mod)-1))
    col_data <- data.frame(col_data, covar_mod[, -1])
    rownames(col_data) <- colnames(counts_mat)
    design_formula <- as.formula(paste0("~", paste(colnames(col_data), collapse="+")))
  }

  dds <- DESeqDataSetFromMatrix(countData=counts_mat, colData=col_data, design=design_formula)
  dds <- DESeq(dds)
  de_res <- results(dds, name="Group_1_vs_0")

  de_called <- rownames(de_res)[which(de_res$pvalue < alpha.unadj)]
  de_called_fdr <- rownames(de_res)[which(de_res$padj < alpha.fdr)]
  return(list(unadj=de_called, fdr=de_called_fdr, de_res=de_res, design=design_formula))
}
