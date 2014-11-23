# Utility functions for high-dimensional cox models
# 
# Author: mg14
###############################################################################



#' Test for interactions
#' A coxph() model is used to systematically test for interactions. Returned are both a Wald test and a likelihood ratio test, 
#' as well as the estimated coefficient and an indicator if the fitting caused a warning(). 
#' @param data A data.frame() with covariates
#' @param survival A Surv()ival object
#' @param pairs A list() with interactions to be considered
#' @param whichMain A vector() of main effects to include. The main effects of the pair of interests are always included/
#' @param minObs The minimal number of observations for main effects other than the interaction. 
#' @param mc.cores The number of cores to use.
#' @return A data.frame() with columns pWald, pLR, coef and warn.
#' 
#' @author mg14
#' @export
TestInteractionsCPH <- function(data, survival, pairs, whichMain = colnames(data), mc.cores=1, minObs = 5){
	whichMain <- names(which(colSums(data!=0,na.rm=TRUE)[whichMain] >= minObs))
	data <- data + 0 ## Make numeric
	result <- mclapply(pairs, function(s){
				a <- s[1]
				b <- s[2]
				if(length(a)==0 | length(b)==0 | sum(data[,a] & data[,b], na.rm=TRUE) < minObs)
					return(c(1,1,NA,1))
				main <- "." #ifelse(!includeAllMain, paste(a,b, sep="+"), ".") 
				interaction <- paste(a,b, sep=":")
				ix <- union(whichMain, c(a, b))
				X <- data[,ix]
				withCallingHandlers(
						tryCatch({
									warn <- 0
									pWald <- NA
									coef <- NA
									pLR <- NA
									fit1 <- coxph(as.formula(paste("survival ~ ", interaction , " + ",main,sep="" )), data=X)
									coef <- coef(fit1)[interaction]
									idx <- which(names(coef(fit1))==interaction)
									pWald <- pchisq(coef(fit1)[interaction]^2/diag(fit1$var)[idx], 1, lower.tail=FALSE)
									fit2 <- coxph(as.formula(paste("survival ~ ", main, sep="" )), data=X)
									pLR <- pchisq(-2 *(fit2$loglik[2] - fit1$loglik[2]), 1, lower.tail=FALSE)
									c(pWald, pLR, coef, warn)
								}, 
								error = function(e) rep(NA,4)),
						warning = function(w) {
							warn <<- 1
							invokeRestart("muffleWarning")
						}
				)
			}, mc.cores=mc.cores)
	result <- as.data.frame(t(as.data.frame(result)))
	rownames(result) <-  sapply(pairs,paste,collapse=":")
	names(result) <- c("pWald","pLR","coef","warn")
	return(result)
}


#' Pairwise interactions.
#' This function creates all pairwise interaction (product) terms of two data.frames()
#' @param X 
#' @param Y 
#' @return A data.frame() of dimensions nrow(X) by (ncol(X) x ncol(Y))
#' 
#' @author mg14
#' @export
MakeInteractions <- function(X,Y){
	Z <- do.call(cbind, lapply(X, `*`, Y))
	colnames(Z) <- apply(expand.grid(colnames(Y), colnames(X)),1,paste, collapse=":")
	return(Z)
}

#' Convert factor to integer.
#' @param F A factor
#' @return A data.frame() with columns corresponding to levels() in the factor. 
#' 
#' @author mg14
#' @export
MakeInteger <- function(F){
	res <- as.data.frame(lapply(levels(F), `==`, F))
	colnames(res) <- levels(F)
	res + 0
}

#' Split a data.frame by a time-dependent variable for use in a coxph model
#' @param dataFrame A data frame
#' @param timeEvent The time of the time-deependent covariate (assumed to be a simple event, such as a transplant).
#' @param timeStop Start time for the survival object. Default: rep(0, nrow(dataFrame)).
#' @param timeStart Stop time of the survival object. Default: dataFrame$time.
#' @param status Status of the survival object (0=alive, 1=dead). Default: dataFrame$status.
#' @return A data.frame with extra extra rows for observations after the event and extra columns time1, time2 (start and stop times), event for the event.
#' 
#' @author mg14
#' @export
MakeTimeDependent <- function(dataFrame, timeEvent, timeStop = dataFrame$time, timeStart = rep(0, nrow(dataFrame)), status=dataFrame$status){
	w <- which(timeEvent < timeStop & timeEvent > timeStart)
	index <- c(1:nrow(dataFrame), w) 
	d <- dataFrame[index,]
	d$index <- index
	d$time1 <- c(timeStart, timeEvent[w])
	d$time2 <- c(pmin(timeStop, timeEvent, na.rm=TRUE), timeStop[w])
	d$event <- c(rep(0,nrow(dataFrame)), rep(1, length(w)))
	e <- c(status, status[w])
	e[w] <- 0
	d$status <- e
	return(d)
}