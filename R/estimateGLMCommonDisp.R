estimateGLMCommonDisp <- function(y, ...)
UseMethod("estimateGLMCommonDisp")

estimateGLMCommonDisp.DGEList <- function(y, design=NULL, method="CoxReid", subset=10000,
                                          verbose=FALSE, lambda_reg=0, alpha_reg=0, ...)
{
#	Check y
	y <- validDGEList(y)
	AveLogCPM <- aveLogCPM(y, dispersion=0.05, lambda_reg=lambda_reg, alpha_reg=alpha_reg)

	disp <- estimateGLMCommonDisp(y=y$counts, design=design, offset=getOffset(y),
	                              method=method, subset=subset, AveLogCPM=AveLogCPM,
	                              verbose=verbose, weights=y$weights, lambda_reg=lambda_reg,
	                              alpha_reg=alpha_reg, ...)

	y$common.dispersion <- disp
	y$AveLogCPM <- aveLogCPM(y, dispersion=disp, lambda_reg=lambda_reg, alpha_reg=alpha_reg)
	y
}

estimateGLMCommonDisp.default <- function(y, design=NULL, offset=NULL, method="CoxReid",
                                          subset=10000, AveLogCPM=NULL, verbose=FALSE,
                                          weights=NULL, lambda_reg=0, alpha_reg=0, ...)
{
#	Check y
	y <- as.matrix(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}
	if(ncol(design) >= ncol(y)) {
		warning("No residual df: setting dispersion to NA")
		return(NA_real_)
	}

#	Check method
	method <- match.arg(method, c("CoxReid","Pearson","deviance"))
	if(!method == "CoxReid" && !is.null(weights)) warning("weights only supported by CoxReid method")

#	Check offset
	if(is.null(offset)) offset <- log(colSums(y))

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y, offset=offset, weights=weights, lambda_reg=lambda_reg, alpha_reg=alpha_reg)

#	Call lower-level function
	disp <- switch(method,
		CoxReid=dispCoxReid(y, design=design, offset=offset, subset=subset, AveLogCPM=AveLogCPM, weights=weights, lambda_reg=lambda_reg, alpha_reg=alpha_reg, ...),
		Pearson=dispPearson(y, design=design, offset=offset, subset=subset, AveLogCPM=AveLogCPM, ...),
		deviance=dispDeviance(y, design=design, offset=offset, subset=subset, AveLogCPM=AveLogCPM, ...)
	)

	if(verbose) cat("Disp =",round(disp,5),", BCV =",round(sqrt(disp),4),"\n")
	disp
}

