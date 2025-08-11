## reComBat-seq

reComBat-seq is a batch effect adjustment tool designed for bulk RNA-seq count data with underdetermined experimental designs. Building upon the negative binomial regression framework used in ComBat-seq[2], reComBat-seq incorporates Elastic Net regularization to address convergence issues when handling confounded data. It takes raw **untransformed, raw** count matrices as input and requires a known batch variable.

By applying regularized negative binomial regression, reComBat-seq models batch effects while enabling stable parameter estimation with rank-deficient design matrices. The adjusted data produced preserves the integer nature of counts, ensuring compatibility with widely used differential expression tools such as edgeR[1] and DESeq2. 

This formulation extends the standard negative binomial regression model by incorporating regularization penalties, with `NLL`representing the negative log-likelihood function.
```math
L(\beta) = NLL(\beta) + \lambda \bigg( \alpha \| \beta \|_1 + \bigg( \frac{1-\alpha}{2} \bigg) \| \beta \|_2^2 \bigg)
```

## Installation

reComBat-seq is available for installation via GitHub.

```r
# install.packages("devtools")
devtools::install_github("jas-st/reComBat-seq")
```


## Usage and Tutorial

An example use case has been included in the `tutorial` folder. It analyses an Esophagus tissue dataset, containing cancer and healthy cells. 
A raw count matrix from RNA-Seq studies has to be provided, without any normalization or transformation, as well as a vector for batch separation). The `group` parameter specifies additional biological covariates, in this case the disease label of the cells.

```r
recombatseq_df <- reComBat_seq(counts_df_reduced, batch = batches, group = group,
                             lambda_reg=0.8, alpha_reg=0.3)
```

| Raw Data      | Corrected Data |
| ------------- | -------------  |
| <img src="https://github.com/jas-st/reComBat-seq/blob/main/tutorial/PCA_raw.png" width="500">  | <img src="https://github.com/jas-st/reComBat-seq/blob/main/tutorial/PCA_corrected.png" width="500">   |
  
For multiple biological variables the `covar_mod` parameter can be used. In this scenario the confounded matrix from the tutorial will be used.

```r
recombatseq_df_conf <- reComBat_seq(counts_df_reduced, batch = batches, covar_mod = covmat_model,
                                   lambda_reg=0.8, alpha_reg=0.3)
```

| Corrected Data      | Corrected Data (Singular matrix) |
| ------------- | -------------  |
| <img src="https://github.com/jas-st/reComBat-seq/blob/main/tutorial/PCA_corrected.png" width="500">  | <img src="https://github.com/jas-st/reComBat-seq/blob/main/tutorial/PCA_corrected_confounded.png" width="500">   |

## Arguments

Since reComBat-seq follows the same structure as ComBat-seq, the parameters used are also the same. The regularization parameters added are the following:

  - `lambda_reg` - controls the strength of the regularization, $\lambda$ in the above equation
  - `alpha_reg` - controls the elastic net tuning, $\alpha$ in the above equation

## References
1. Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2025). “edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.” Nucleic Acids Research, 53(2), gkaf018. doi:10.1093/nar/gkaf018.
2. Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078
