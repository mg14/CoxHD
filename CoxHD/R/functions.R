# Utility functions for high-dimensional cox models
# 
# Author: mg14
###############################################################################


##### Simulation study for ecoxph
SimSurv <- function(risk) {
	n = length(risk)
	death = 1/log(2) * log(rexp(n, exp(risk)) * log(2) +1)
	cens = rbinom(n,1,0.5)
	surv = Surv(time = pmax(0,death * pmax(runif(n,0.5,1),1-cens)), event=1-cens)
	return(surv)
}

#' Non-parametric survival simulations
#' @param risk A vector of the relative risk/
#' @param surv A survival object
#' @param H0 The baseline hazed function. By default this is given by basehaz(coxph(surv ~ 1)).
#' @return A Surv()ival object
#' 
#' @author mg14
#' @export
SimSurvNonp <- function(risk, surv, H0 = basehaz(coxph(surv ~ 1))) {
	## Simulate deaths times
	FhazInv <- splinefun(c(0,H0$hazard), c(0,H0$time), method="monoH.FC")
	n = length(risk)
	deathTimes = FhazInv(rexp(n, exp(risk))) #predict(hazardDist, rexp(n, exp(risk)))
	
	## Simulate censoring times
	nCens <- sum(surv[,2]==0, na.rm=TRUE)
	Fsurv <- splinefun(exp(-H0$hazard) ~ H0$time, method="monoH.FC")
	x <- sort(na.omit(surv[surv[,2]==0,1])) ## observed (conditioned) censored times
	t <- table(x) / nCens
	y <- t / Fsurv(unique(x)) ## Need to get unconditional distribution
	Fcens <- cumsum(y)
	Fcens <- Fcens/max(Fcens)
	FcensInv <- splinefun(unique(x) ~ Fcens, method="monoH")
	censTimes <- FcensInv(runif(n,0,1))
	
	## Put together
	survOut <-  Surv(time = pmax(0,pmin(deathTimes, censTimes)), event=(deathTimes < censTimes)+0)
	return(survOut)
}


#' Non-parametric data extrapolations based on multiple imputation
#' A fraction percentMissing is set to NA and multiple imputation from the mice package is used to fill these gaps.
#' @param oldData 
#' @param nData 
#' @param percentMissing 
#' @param ... 
#' @return A data.frame with extrapolatied data
#' 
#' @author mg14
#' @export
SimDataNonp <- function(oldData, nData, percentMissing = 0.33, ...){
	require(mice)
	oldData <- oldData[!apply(is.na(oldData), 1, all),]
	for(i in 1: nrow(oldData)){
		while(TRUE){
			naIdx <- sample(ncol(oldData), round(percentMissing * ncol(oldData)))
			if(! all(1:ncol(oldData) %in% naIdx))
				break
		}
		oldData[i,naIdx] <- NA
	}
	#oldData <- as.matrix(oldData[sample(nrow(oldData), nData, replace=TRUE),])
	m <- mice(as.data.frame(oldData), printFlag=FALSE, ...)
	newData <- complete(m, action="long")
	newData <- newData[sample(nrow(newData), nData, replace=nrow(newData)< nData),-2:-1]
	return(newData)
}

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