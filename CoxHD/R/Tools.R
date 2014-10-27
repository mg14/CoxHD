# Utility functions for high-dimensional cox models
# 
# Author: mg14
###############################################################################


GetPairs <- function(names, scope){
	pairs <- strsplit(scope, "\\.")
	lapply(pairs, function(s){
				a <- grep(paste(s[1],"$", sep=""), names, value=TRUE)
				b <- grep(paste(s[2],"$", sep=""), names, value=TRUE)
				c(a,b)
			})
}

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
TestInteractions <- function(data, survival, pairs, whichMain = colnames(data), mc.cores=1, minObs = 5){
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

DensityEstimates <- function(coxRFX, newx = range(coef(coxRFX)), n = 100){
	x <- seq(newx[1], newx[2], length.out=n)
	c <- coef(coxRFX)
	v <- diag(coxRFX$var)
	z <- mapply(function(i,j) dnorm(x, i, j), c, sqrt(v))
	sapply(levels(coxRFX$groups), function(g) rowMeans(z[,coxRFX$groups==g, drop=FALSE]))
}


#### Paralllelized stepwise forward selection
parStep = function(X, surv, criterion = "AIC", max.iter = ncol(X), mc.cores=1, verbose=FALSE){
	select = numeric()
	loglik = numeric()
	penalty = ifelse(criterion == "AIC",1,log(nrow(X))/2)
	while(length(select) < max.iter){
		k = length(select) + 1
		scope = setdiff(1:ncol(X), select)
		if(is.null(scope))
			break
		l = mclapply(scope, function(i) {
					x = as.data.frame(X[,c(select,i)])
					t = try(coxph(surv~., data=x)$loglik[2])
					ifelse(class(t)!="try-error",t,-Inf)
				}, mc.cores=mc.cores)
		logliks = unlist(l)		
		add = scope[which.max(logliks)]
		if(k>1)
			if(all(logliks <= loglik[k-1] + penalty))
				break
		loglik = c(loglik,max(logliks))
		select = c(select,add)
		if(verbose) cat(".")
	}
	if(verbose) cat("\n")
	return(data.frame(select=select,loglik=loglik,AIC = - 2*loglik + 2 * 1:length(loglik), BIC = - 2 * loglik +  1:length(loglik) * log(nrow(X)) ))
}

#### Convert Pi to p-values
pi2p = function(s, pi.thr=0.5){
	q = colMeans(s$Pr)
	p.conc = sapply(q, minD, B=s$bootstrap.samples)
	Lambda = apply(p.conc, 2, function(x) which(x < 0.05)[1] )/(2*s$bootstrap.samples) < pi.thr
	pi = apply(s$Pr,2,max)
	p.val = c(1,minD(mean(pi), B = s$bootstrap.samples))[pi * 2 * s$bootstrap.samples +1]
	return(data.frame(pi=pi, p.val=p.val))
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

MakeTimeDependent <- function(dataFrame, timeTpl, timeSurv = dataFrame$time, time0Surv = rep(0, nrow(dataFrame)), event=dataFrame$event){
	w <- which(timeTpl < timeSurv & timeTpl > time0Surv)
	index <- c(1:nrow(dataFrame), w) 
	d <- dataFrame[index,]
	d$index <- index
	d$time1 <- c(time0Surv, timeTpl[w])
	d$time2 <- c(pmin(timeSurv, timeTpl, na.rm=TRUE), timeSurv[w])
	d$transplant <- c(rep(0,nrow(dataFrame)), rep(1, length(w)))
	e <- c(event, event[w])
	e[w] <- 0
	d$event <- e
	return(d)
}