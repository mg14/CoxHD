# A random effects Cox model
# 
# Author: mg14
###############################################################################


#' 
#' @param X 
#' @param surv 
#' @param tol 
#' @param max.iter 
#' @return NA
#' 
#' @author mg14
#' @noRd
ecoxph <- function(X,surv, tol=1e-3, max.iter=50){
	if(class(X)=="data.frame")
		X = as.matrix(X)
	beta0 = rep(0,ncol(X))
	beta1 = rep(1,ncol(X))
	sigma2 = 1
	iter = 1
	while(max(abs(beta1-beta0))>tol& iter < max.iter){
		fit = coxph(surv ~ ridge(X, theta=1/sigma2, scale=FALSE))
		sigma2 = (1 + sum((fit$coefficients-mean(fit$coefficients))^2))/(ncol(X))	
		beta0 = beta1
		beta1 = fit$coefficients
		#cat(beta1,"\n")
		#cat(sigma,"\n")
		iter = iter+1
	}
	fit$sigma2 = sigma2
	names(fit$coefficients) = colnames(X)
	return(fit)
}

#' Cox proportional hazards model with random effects
#' 
#' This functinon estimates a Cox proportional in which the parameters follow normal distributions. Multiple groups can be defined with different prior mean and variance. 
#' The parameters of the prior distribution are estimated by an empirical Bayes approach.
#' The model extends the shared frailty model implemented in the survival package to the situation where the frailty covariates are not a factor, but can be overlapping groups.  
#' @param X The data matrix (n x p)
#' @param surv The survival object (n x 2)
#' @param groups Optional groups as a factor (p) with l levels. Default = rep(1, n)
#' @param which.mu Indicator which of the groups should have an offset. Default = unique(groups)
#' @param tol The tolerance beyond which to stop
#' @param max.iter The maximal number of iterations
#' @param sigma0 The variance of a si-chisq hyperprior on the variances.
#' @param nu The df of the variance hyperprior. Default = 0, that is no hyperprior.
#' @param penalize.mu Wether to define an N(0,tau) hyperprior on the group means.
#' @param sigma.hat Which estimator to use for the variances. Default MLE. Other possibilities include REML and BLUP. 
#' @param verbose Gives more output.
#' @return A coxph object with a few extra fields: $groups, $X, $surv, $sigma2 (the variances), $mu (the means)
#' 
#' @author mg14
#' @export
CoxRFX <- function(X, surv, groups = rep(1, ncol(X)), which.mu = unique(groups), tol=1e-3, max.iter=50, sigma0 = 0.1, nu = 0,  penalize.mu = FALSE, sigma.hat=c("MLE","REML","df","BLUP"), verbose=FALSE){
	if(class(X)=="data.frame")
		X = as.matrix(X)
	sigma.hat = match.arg(sigma.hat)
	o <- order(groups)
	X <- X[,o]
	groups <- factor(groups[o])
	uniqueGroups <- levels(groups)
	XX <- lapply(uniqueGroups, function(i) X[,groups==i, drop=FALSE])
	names(XX) <- uniqueGroups
	sumX <- sapply(which.mu, function(i) rowSums(XX[[i]]))
	nGroups = length(uniqueGroups)
	sigma2 <- sigma0ld <- rep(ifelse(sigma0>0, sigma0,1), nGroups)
	iter = 1
	mu <- mu0ld <- rep(0, nGroups)
	names(mu) <- uniqueGroups
	beta = rep(1,ncol(X)+length(which.mu))
	beta0ld = rep(0,ncol(X)+length(which.mu))
	sigma2.mu = 42
	if(!is.null(which.mu)) 
		if(!penalize.mu)
			sumTerm <- "sumX" 
		else
			sumTerm <- "ridge(sumX, theta=1/sigma2.mu, scale=FALSE)"
	else sumTerm <- character(0)
	while((max(abs(beta-beta0ld)) > tol | max(abs(mu - mu0ld)) > tol | max(abs(sigma2 - sigma0ld)) > tol) & iter < max.iter){
		beta0ld = beta
		sigma0ld <- sigma2
		mu0ld <- mu
		formula <- formula(paste("surv ~", paste(c(sapply(1:nGroups, function(i) paste("ridge(XX[[",i,"]], theta=1/sigma2[",i,"], scale=FALSE)", sep="")), 
								#ifelse(!is.null(which.mu),"ridge(sumX, theta=1/sigma.mu, scale=FALSE)","")), 
								sumTerm), 
						collapse=" + ")))
		fit <- coxph(formula)
		if(!is.null(which.mu))
			mu[which.mu] <- coef(fit)[-(1:ncol(X))]
		if(verbose) cat("mu", mu, "\n", sep="\t")
		names(fit$df) <- c(uniqueGroups, rep("Offset", length(which.mu)>0))
		if(verbose) cat("df", fit$df,"\n", sep="\t")
		sigma2 = sapply(uniqueGroups, function(i){
					index <- which(groups==i) #& fit$coefficients > beta.thresh
					if(sigma.hat=="BLUP")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ))/(nu + length(index) -1) #+ mean(diag(fit$var)[index])
					else if(sigma.hat=="df")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ))/(nu + fit$df[i]) #+ mean(diag(fit$var)[index]) ## 
					else if(sigma.hat == "MLE")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ) + sum(diag(solve(solve(fit$var)[index,index]))))/(nu + length(index))
					else if(sigma.hat == "REML")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ) + sum(diag(fit$var)[index]))/(nu + length(index))
				})
		if(verbose) {
			cat("sigma2", sigma2, "\n", sep="\t")
			cat("loglik:", fit$loglik - c(0,fit$penalty[2] + 1/2 * sum(log(sigma2[groups]))),"\n", sep="\t")
		}
		if(penalize.mu){
			if(sigma.hat=="BLUP")
				sigma2.mu = (sigma0 * nu + sum((mu-0)^2)) / (nu + length(mu))
			else if(sigma.hat=="df")
				sigma2.mu = (sigma0 * nu + sum((mu-0)^2)) / (nu + fit$df["Offset"])
			else if(sigma.hat == "MLE")
				sigma2.mu = (nu * sigma0 + sum((mu-0)^2 ) + sum(diag(solve(solve(fit$var)[-(1:ncol(X)),-(1:ncol(X))]))))/(nu + length(mu))
			else if(sigma.hat == "REML")
				sigma2.mu = (nu * sigma0 + sum((mu-0)^2 ) + sum(diag(fit$var)[-(1:ncol(X))]))/(nu + length(mu))
		}
		
		#cat(sigma.mu,"\n")
		beta = fit$coefficients
		
		#beta1[beta1 < beta.thresh] <- 0
		#c = lapply(unique(groups), function(i) beta1[groups==i])
		#for(i in 1:nGroups)
		#	XX[[i]][,c[[i]]==0] <- 0
		#cat(beta1,"\n")
		#cat(sigma,"\n")
		#cat(max(abs(beta - beta0ld)), max(abs(mu - mu0ld)), max(abs(sigma2 - sigma0ld)), "\n", sep="\t")
		iter = iter+1
	}
	if(iter == max.iter)
		warning("Did not converge after", max.iter, "iterations.")
	fit$iter[1] <- iter
	fit$sigma2 = sigma0ld
	names(fit$sigma2) <- uniqueGroups
	fit$sigma2.mu = sigma2.mu
	fit$mu = mu
	#fit$sumX = sumX
	fit$X = X[,order(o)]
	fit$surv = surv
	fit$groups = groups[order(o)]
	var = fit$var
	var2 = fit$var2
	fit$var = var[1:ncol(X),1:ncol(X)][order(o),order(o)]
	fit$var2 = var2[1:ncol(X),1:ncol(X)][order(o),order(o)]
	fit$mu.var = var[-(1:ncol(X)),-(1:ncol(X))]
	fit$mu.var2 = var2[-(1:ncol(X)),-(1:ncol(X))]
	fit$means = fit$means[1:ncol(X)][order(o)]
	#fit$delta = sapply(unique(groups), function(i) mean(groups==i & fit$coefficients < beta.thresh))
	fit$coefficients <- fit$coefficients[1:ncol(X)][order(o)] + mu[fit$groups]
	names(fit$coefficients) = colnames(X)[order(o)]
	fit$terms <- fit$terms[1:length(uniqueGroups)]
	fit$penalized.loglik <- fit$loglik[2] - fit$penalty[2] - 1/2 * sum(log(fit$sigma2[groups]))
	class(fit) <- c("CoxRFX", class(fit))
	return(fit)
}

#' Partial risk components
#' @param fit The CoxRFX fit
#' @param newX New data, defaults to fit$X
#' @param groups The groups, defaults to fit$groups
#' @return A matrix with the risk per group
#' 
#' @author mg14
#' @export
PartialRisk <- function(fit, newX=fit$X, groups=fit$groups) {
	sapply(levels(groups), function(x) {
				ix <- groups == x
				as.matrix(newX[,ix, drop=FALSE]) %*% coef(fit)[ix] 
			})
}

#' Variance (confidence intervals) of partial risk components
#' @param fit The CoxRFX fit
#' @param newX New data, defaults to fit$X
#' @param groups The groups, defaults to fit$groups
#' @return A matrix with the confidence interval (prediction variance) for each risk component
#' 
#' @author mg14
#' @export
PartialRiskVar <- function(fit, newX=fit$X, groups=fit$groups) {
	newX <- newX - rep(colMeans(newX), each=nrow(newX))
	sapply(levels(groups), function(x) {
				ix <- groups == x
				rowSums((as.matrix(newX[,ix, drop=FALSE]) %*% fit$var[ix,ix]) * as.matrix(newX[,ix, drop=FALSE]))
			})
}

#' Variance components
#' @param fit The CoxRFX fit
#' @param newX New data, defaults to fit$X
#' @param groups The groups, defaults to fit$groups
#' @param type Take either the diagonal of the covariance matrix (default), or the rowSums. The latter guaranties that the
#' components sum up to the actual variance, but could be negative in the case of collinearity. The two are equivalent 
#' @return A vector containing the variance components
#' 
#' @author mg14
#' @export
VarianceComponents <- function(fit, newX = fit$X, groups = fit$groups, type = c("diag","rowSums")){
	risk <- PartialRisk(fit = fit, newX = newX, groups = groups)
	type <- match.arg(type)
	#residual <- predict(fit, se.fit=TRUE)$se.fit^2
	newX <- as.matrix(newX - rep(colMeans(newX), each=nrow(newX)))
	residual <- rowSums((newX %*% fit$var) * newX)
	
	c <- cov(risk, use="complete")
	if(type=="diag")
		x <- diag(c)
	else
		x <- rowSums(c)
	return(c(x, residual=mean(residual)))
}

#' Plot variance components
#' @param fit The CoxRFX fit
#' @param col The colors for each component
#' @param groups the groups to be used, if different from the fitted ones.
#' @param type The type of variance compnents: Either 'rowSums' (default) or 'diag'. Rowsums sum up to the actual variance of the linear predictor, but can be
#' negative. Plotting just the 'diag'onal elements guarantees positive components  
#' @return NULL
#' 
#' @author mg14
#' @export
PlotVarianceComponents <- function(fit, col=1:nlevels(fit$groups), groups = fit$groups, type="rowSums", conf.int=TRUE, digits=2) {
	if(is.null(names(col)))
		names(col) <- levels(groups)
	v <- VarianceComponents(fit, groups=groups, type=type)
	o <- order(v[levels(groups)], decreasing=TRUE)
	v <- v[o]
	vp <- v[v>0]
	vn <- v[v<=0]
	if(conf.int){
		r <- paste( "+/-",round(colMeans(PartialRiskVar(fit, groups=groups))[o], digits))
		rp <- r[v>0]
		rn <- r[v<=0]
	}else{
		rp <- rn <- NULL
	}
	pie(vp, col=col[names(vp)], border=NA, labels=paste(names(vp), " (", round(vp, digits),rp,")",sep=""), radius = sum(v))
	
	if(length(vn)>0){
		par(new=T)
		pie(c(abs(vn), sum(vp)-sum(vn)), col=c(col[names(vn)],NA), border=NA, labels=paste(names(vn), " (", round(vn, digits),rn,")", sep=""), new=FALSE, density=c(rep(36, length(vn) ),NA), radius = sum(v))
	}
}

VarianceComponentsCV <- function(fit, which.coef = grep(":", colnames(fit$X), invert = TRUE, value=TRUE), type=c("diag","colSums"), method="simple"){
	if(method=="simple"){
		type="diag"
		risk <- PartialRisk(fit = fit)
		X <- fit$X
		risk0 <- X * rep(fit$coef, each=nrow(X))
		X <- X - rep(colMeans(X), each=nrow(X))
		#residual <- predict(fit, se.fit=TRUE)$se.fit^2
		V <- fit$var
		colnames(V) <- rownames(V) <- colnames(fit$X)
		sapply(which.coef, function(i){
					w <- grep(i, colnames(X), value = TRUE)
					r <- risk
					for(ww in w)
						r[,fit$groups[ww]] <- r[,fit$groups[ww]] - risk0[,ww]
					idx <- !colnames(X) %in% w
					residual <- rowSums((as.matrix(X[,idx]) %*% V[idx,idx]) * as.matrix(X[,idx]))
					c <- cov(r, use="complete")
					if(type=="diag")
						x <- diag(c)
					else
						x <- rowSums(c)
					return(c(x, residual=mean(residual)))
				})
	}else{
		sapply(which.coef, function(i){
					w <- grep(i, colnames(fit$X), value = TRUE, invert=TRUE)
					fit <- CoxRFX(fit$X[,w], fit$surv, fit$groups[w], sigma0=0.1, nu=0)
					VarianceComponents(fit)
	})}
}


#' Compute concordance for risk components
#' @param fit The CoxRFX fit
#' @param newX New data, defaults to fit$X
#' @param groups The groups, defaults to fit$groups
#' @param newSurv The survival object, defaults to fit$surv
#' @return A vector with the concordance of each component
#' 
#' @author mg14
#' @export
PartialC <- function(fit, newX = fit$X, newSurv=fit$surv, groups=fit$groups){
	#require(Hmisc)
	c <- sapply(levels(groups), function(x) {
				ix <- groups == x
				r <- as.matrix(newX[,ix, drop=FALSE]) %*% coef(fit)[ix]
				#rcorr.cens(-r, newSurv)[1]
				c <- survConcordance(newSurv~r)
				c(c$concordance, c$std.err)
			})
	colnames(c) <- levels(groups)
	return(c)
}
	
#' Predict risk with missing data
#' @param fit A CoxRFX fit
#' @param newX 
#' @param var Which variance estimate to use either var = $H^{-1}$, or var2 = $H^{-1}I H^{-1}$. The former is the more conservative choice, the latter seems more accurate, but may underestimate the variance.
#' @return A data.frame with columns Expected and Variance 
#' 
#' @author mg14
#' @export
PredictRiskMissing <- function(fit, newX=fit$X, var = c("var","var2")){
	var <- match.arg(var)
	Sigma <- cov(fit$X)
	mu <- colMeans(fit$X)
	beta <- fit$coefficients
	
	.predict <- function(newX, beta, Sigma, mu){
		missing <- is.na(newX)
		expectedX <- newX
		if(any(missing)){
			s <- Sigma[missing, !missing] %*% MASS::ginv(Sigma[!missing, !missing])
			expectedX[missing] <- mu[missing] + s %*% (newX[!missing] - mu[!missing])
			varianceRisk <- beta[missing] %*% (Sigma[missing,missing] - s %*%  Sigma[!missing, missing] ) %*% beta[missing]
		}else{
			varianceRisk <- 0
		}
		expectedRisk <- expectedX %*% beta
		e <- expectedX - mu
		varianceRisk <- e %*% fit[[var]][1:length(fit$coefficients),1:length(fit$coefficients)] %*% e + varianceRisk
		return(c(expectedRisk, varianceRisk))
	}
	predictions <- t(apply(newX, 1, .predict, beta, Sigma, mu))
	colnames(predictions) <- c("Expected","Variance")
	return(predictions)
}

#' Covariance-based imputation of missing variables
#' 
#' This functinon imputes missing variables based on the covariance of a previous data set.
#' This can be useful if a CoxRFX model (or any other) has been fit on a large data set with
#' correlated variables and one wants to use this model on a different data sets with fewer 
#' or randomly missing covariates. 
#' @param X orignial data set
#' @param newX The data.frame of covariates
#' @return A data.frame of dim(newX) with imputed variables
#' 
#' @author mg14
#' @export
ImputeXMissing <- function(X, newX){
	Sigma <- cov(X)
	mu <- colMeans(X)
	l <- ncol(X)
	
	.impute <- function(newX, Sigma, mu){
		missing <- is.na(newX)
		expectedX <- newX
		varianceX <- rep(0,l)
		if(any(missing)){
			s <- Sigma[missing, !missing] %*% MASS::ginv(Sigma[!missing, !missing])
			varianceX[missing] <- diag(s)
			expectedX[missing] <- mu[missing] + s %*% (newX[!missing] - mu[!missing])
		}
		#return(cbind(expectedX, varianceX))
		expectedX
	}
	imputations <- t(apply(newX, 1, .impute, Sigma, mu))
	colnames(imputations) <- colnames(newX) #c("Expected","Variance")
	return(imputations)
}

#' Standardize the magnitude of covariates
#' @param X 
#' @return data.frame of dim(X)
#' 
#' @author mg14
#' @export
StandardizeMagnitude <- function(X){
	.scale <- function(x) {if(max(x, na.rm=TRUE)==0) 1 else 10^floor(log10(max(x, na.rm=TRUE)))}
	scale <- sapply(X, .scale)
	Y <- X * rep(1/scale, each=nrow(X))
	n <- as.character(scale)
	i <- n!="1"
	names(Y)[i] <- paste(names(X)[i], n[i], sep="_")
	return(Y)
}

#' Plot a CoxRFX model
#' @param fit 
#' @param col 
#' @param order 
#' @param xlim 
#' @param xlab 
#' @param ... 
#' @return NULL
#' 
#' @author mg14
#' @export
plot.CoxRFX <- function(fit, col=c(brewer.pal(9,"Set1"), brewer.pal(8,"Dark2")), order = 1:nlevels(fit$groups), xlim=range(coef(fit)), xlab="Coefficient",...){
	plot(NA,NA, xlim=xlim, ylim=range(1,1+nlevels(fit$groups)), yaxt="n", ylab="",xlab=xlab, ...)
	axis(side=2, at=1:nlevels(fit$groups), labels = levels(fit$groups)[order], las=2)
	i <- 1
	for(l in levels(fit$groups)[order]){
		m <- fit$mu[l]
		s <- fit$sigma2[l]
		x <- seq(m-3*sqrt(s),m+3*sqrt(s), l=100)
		y <- i+dnorm(x, m, sqrt(s))/dnorm(m, m, sqrt(s))*.8
		polygon(x,y, col=paste(col[order][i],"44", sep=""), border=NA)
		lines(x, y, col=col[order][i])
		points(fit$coef[fit$groups==l],rep(i, sum(fit$groups==l)), col=col[order][i], pch=16, cex=.5)
		lines(rep(m,2),c(0,.8) +i, col=col[order][i])
		i <- i+1
	}
}
