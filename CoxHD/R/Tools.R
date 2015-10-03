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
	d$time2 <- c(timeStop, timeStop[w])
	d$time2[w] <- timeEvent[w]
	d$event <- c(rep(0,nrow(dataFrame)), rep(1, length(w)))
	e <- c(status, status[w])
	e[w] <- 0
	d$status <- e
	return(d)
}


#' Post-hoc competing risk adjustment for two Kaplan-Meier curves.
#' 
#' The function uses a monotonous spline fit to interpolate between the steps of the fit.
#' @param fit1 
#' @param fit2 
#' @return a survfit object with updated incidence and confidence intervals
#' 
#' @author mg14
crAdjust <- function(fit1, fit2){
	int2 <- splinefun(fit2$time, fit2$surv,  method="monoH.FC")
	fit1$surv <- cumsum(c(1,diff(fit1$surv)) * int2(fit1$time))
	fit1$upper <- cumsum(c(1,diff(fit1$upper)) * int2(fit1$time))
	fit1$lower <- cumsum(c(1,diff(fit1$lower)) * int2(fit1$time))
}

#' Extract and predict survival status
#' 
#' This function extracts the survival status at a given time t. For censored patients at time t_cens < t, the 
#' predicted status is S[t]/S[t_cens] if censored="conditional" and NA otherwise.
#' @param surv As Surv() object. Currently only right-censoring supported.
#' @param time The time at which to extract the survival status.
#' @param censored A character vector, either "conditional" or "na".
#' @return A vector with the survival status.
#' 
#' @author mg14
#' @export
survStatus <- function(surv, time, censored=c("conditional","na")){
	censored <- match.arg(censored)
	o <- ifelse(ncol(surv)==2, 0, 1)
	if(o==1) warning("Only right-censored data are supported.")
	pSurv <- function(time, surv){
		s <- summary(survfit(surv ~ 1))
		w <- rowSums(matrix(rep(s$time, each=length(time)), nrow=length(time))<=time)
		c(1,s$surv)[w+1]
	}
	absSurv <- pSurv(time, surv)
	status <- surv[,1+o] >= time
	w <- which(surv[,2+o]==0 & surv[,1+o] < time)
	if(length(w) > 0){
		if(censored=="na") 
			condSurv <- NA
		else 
			condSurv <- absSurv/pSurv(surv[w,1+o], surv)
		status[w] <- condSurv
	}
	return(status)
}

#' Evaluation of absolute prediction errors. 
#' 
#' Given predictions p in [0,1], the function evaluates the mean absolute prediction error sum(|p-status|)/n, the mean squared error (Brier score)
#' sum((p-status)^2)/n, the log2 entropy sum(p * log2(status/p))/n, and a Bayes misclassification rate sum((p>.5)*(1-status))/n. If the 
#' status is extrapolated with probability in [0,1], weighted averages are reported.
#' @param x A vector with the predictions in [0,1]
#' @param surv The Surv() object of the observed survival.
#' @param time The time at which the predictions are to be evaluated.
#' @param censored A character vector indicating the strategy to handle cases censored before time.
#' @return A vector of length
#' 
#' @author mg14
#' @export
ape <- function(x, surv, time, censored="conditional"){
	status <- survStatus(surv, time, censored=censored)
	err <- c(abs=mean(status*(1-x) + (1-status)*x, na.rm=TRUE),
			brier=mean(status*(1-x)^2 + (1-status)*x^2, na.rm=TRUE),
			log2=mean((log2((status/x)^status)+log2(((1-status)/(1-x))^(1-status))), na.rm=TRUE),
			bayes=mean((x>.5) != (status>.5), na.rm=TRUE))
	return(err)
}