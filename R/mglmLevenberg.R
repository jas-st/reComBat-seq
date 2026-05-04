# Adapted from edgeR 4.0.1 (mglmLevenberg.R)
# Modified by Zhasmina Stoyanova, 2026
# Changes: added lambda_reg, alpha_reg, and num_threads parameters for elastic net
# regularization and OpenMP parallelization, passed through to .cxx_fit_levenberg.
# Original authors: edgeR team (see edgeR package).
mglmLevenberg <- function(y, design, dispersion=0, offset=0, weights=NULL,
                          coef.start=NULL, start.method="null", maxit=200, tol=1e-06,
                          lambda_reg=0, alpha_reg=0, num_threads=0)
#	Fit genewise negative binomial glms with log-link
#	using Levenberg damping to ensure convergence

#	R version by Gordon Smyth and Yunshun Chen
#	C++ version by Aaron Lun
#	Created 3 March 2011.  Last modified 21 June 2017
{

# Save current state to restore it later
	if(num_threads) {
	    old_blas_threads <- RhpcBLASctl::blas_get_num_procs()
	    
	    # FORCE BLAS TO 1 THREAD
	    # This prevents the "OpenBLAS loop failed" error on compute nodes
	    RhpcBLASctl::blas_set_num_threads(1)
	    
	    # Ensure threads are reset even if the C++ code throws an error or user hits 'Esc'
	    on.exit({
	      RhpcBLASctl::blas_set_num_threads(old_blas_threads)
	    }, add = TRUE)
	    
	    message(
			sprintf("Parallel check: Throttling BLAS (%d -> 1) for %d OpenMP threads", old_blas_threads, num_threads)
		)
		omp_set_num_threads(num_threads)
	}

      
#	Check arguments
	y <- as.matrix(y)
	if(!is.numeric(y)) stop("y is non-numeric")
	nlibs <- ncol(y)
	ngenes <- nrow(y)
	if(nlibs==0 || ngenes==0) stop("no data")

#	Checks for negative, NA or non-finite values in the count matrix.
	.isAllZero(y)

#	Checking the design matrix
	design <- as.matrix(design)
	if (!is.double(design)) storage.mode(design) <- "double"
	if (!all(is.finite(design))) stop("all entries of design matrix must be finite and non-missing")

#	Checking dispersions, offsets and weights
	offset <- .compressOffsets(y, offset=offset)
    dispersion <- .compressDispersions(y, dispersion)
	weights <- .compressWeights(y, weights)

#	Initializing values for the coefficients at reasonable best guess with linear models.
	if(is.null(coef.start)) {
		start.method <- match.arg(start.method, c("null","y"))
		beta <- .Call(.cxx_get_levenberg_start, y, offset, dispersion, weights, design, start.method=="null")
	} else {
		beta <- as.matrix(coef.start)
	}

# 	Checking arguments and calling the C++ method.
	if (!is.double(beta)) storage.mode(beta) <- "double"
	output <- .Call(.cxx_fit_levenberg, y, offset, dispersion,
	                weights, design, beta, tol, maxit, lambda_reg, alpha_reg, num_threads)

#	Naming the output and returning it.
	names(output) <- c("coefficients", "fitted.values", "deviance", "iter", "failed")
	colnames(output$coefficients) <- colnames(design)
	rownames(output$coefficients) <- rownames(y)
	dimnames(output$fitted.values) <- dimnames(y)
	output
}
