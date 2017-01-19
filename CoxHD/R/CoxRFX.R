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
#' This function estimates a Cox proportional in which the parameters follow normal distributions as discussed by Therneau et al. (2003). 
#' Multiple groups can be defined with different prior mean and variance. 
#' The variances of the joint distributions are efficiently estimated by an EM-type algorithm.
#' @param Z The data matrix of random effects (n x p)
#' @param surv The survival object (n x 2)
#' @param groups Optional groups as a factor (p) with l levels. Default = rep(1, n)
#' @param which.mu Indicator which of the groups should have an offset. Default = unique(groups)
#' @param tol The tolerance beyond which to stop
#' @param max.iter The maximal number of iterations
#' @param sigma0 The variance of a si-chisq hyperprior on the variances.
#' @param nu The df of the variance hyperprior. Default = 0, that is no hyperprior.
#' @param penalize.mu Wether to define an N(0,tau) hyperprior on the group means.
#' @param sigma.hat Which estimator to use for the variances. Default df, other possibilities include MLE, REML and BLUP, see details.
#' @param verbose Gives more output.
#' @details The values of the means mu_g are estimated using the rowSums of Z (within in group) as auxillary variables. 
#' 
#' Different estimators exist for the variances sigma2_g: The default is "df", as used by Perperoglou (2014) and introduced by Schall (1991). In the M-step of the algorithm, this uses sigma^2_g = beta_g beta_g^T/df_g, where the degrees 
#' of freedom df_g = tr H_{gg} are the trace of the Hessian matrix over the elements of group g. Alternatives are MLE, REML, and BLUP, as defined by Therneau et al. (2003). 
#' Simulations indicate that the 'df' method is most accurate.
#' 
#' The model is equivalent to coxme(surv ~ (Z1|1) + rowSums(Z1) + (Z2|1) + rowSums(Z2) + ...); the coxme routine numerically optimises the integrated partial likelihood, which may
#' be more accurate, but is computationally expensive.
#' 
#' @references Terry M Therneau, Patricia M Grambsch & V. Shane Pankratz (2003) Penalized Survival Models and Frailty, Journal of Computational and Graphical Statistics, 12:1, 156-175, http://dx.doi.org/10.1198/1061860031365
#' 
#' A. Perperoglou (2014). Cox models with dynamic ridge penalties on time-varying effects of the covariates. Stat Med, 33:170-80. http://dx.doi.org/10.1002/sim.5921
#' 
#' R. Schall (1991). Estimation in generalized linear models with random effects. Biometrika, 78:719-727. http://dx.doi.org/10.1093/biomet/78.4.719

#' @return A coxph object with a few extra fields: $groups, $Z, $surv, $sigma2 (the variances), $mu (the means), $Hinv (the inverse Hessian of the penalised likelihood), $V = Hinv I Hinv, the covariance of all coefficients and means, 
#' $C the map between centred (beta', mu) to beta. 
#' 
#' @author mg14
#' @export
#' @example inst/example/CoxRFX-example.R
CoxRFX <- function(Z, surv, groups = rep(1, ncol(Z)), which.mu = unique(groups), tol=1e-3, max.iter=50, sigma0 = 0.1, nu = 0,  penalize.mu = FALSE, sigma.hat=c("df","MLE","REML","BLUP"), verbose=FALSE){
	if(class(Z)=="data.frame"){
		Z = as.matrix(Z)
		Z.df <- TRUE
	}else
		Z.df <- FALSE
	if(is.null(colnames(Z)))
		colnames(Z) <- make.names(1:ncol(Z))
	sigma.hat = match.arg(sigma.hat)
	o <- order(groups)
	Z <- Z[,o]
	groups <- factor(groups[o])
	uniqueGroups <- levels(groups)
	ZZ <- lapply(uniqueGroups, function(i) Z[,groups==i, drop=FALSE])
	names(ZZ) <- uniqueGroups
	sumZ <- sapply(which.mu, function(i) rowSums(ZZ[[i]]))
	nGroups = length(uniqueGroups)
	sigma2 <- sigma0ld <- rep(ifelse(sigma0>0, sigma0,1), nGroups)
	iter = 1
	mu <- mu0ld <- rep(0, nGroups)
	names(mu) <- uniqueGroups
	beta = rep(1,ncol(Z)+length(which.mu))
	beta0ld = rep(0,ncol(Z)+length(which.mu))
	sigma2.mu = 42
	if(!is.null(which.mu)) 
		if(!penalize.mu)
			sumTerm <- "sumZ" 
		else
			sumTerm <- "ridge(sumZ, theta=1/sigma2.mu, scale=FALSE)"
	else sumTerm <- character(0)
	while((max(abs(beta-beta0ld)) > tol | max(abs(mu - mu0ld)) > tol | max(abs(sigma2 - sigma0ld)) > tol) & iter < max.iter){
		beta0ld = beta
		sigma0ld <- sigma2
		mu0ld <- mu
		formula <- formula(paste("surv ~", paste(c(sapply(1:nGroups, function(i) paste("ridge(ZZ[[",i,"]], theta=1/sigma2[",i,"], scale=FALSE)", sep="")), 
								#ifelse(!is.null(which.mu),"ridge(sumZ, theta=1/sigma.mu, scale=FALSE)","")), 
								sumTerm), 
						collapse=" + ")))
		fit <- coxph(formula)
		if(any(is.na(coef(fit)))){
			warning(paste("NA during estimation (iter: ", iter, ", coef: ", paste(which(is.na(coef(fit)[order(o)])), sep=","), ")", sep=""))
			break
		}
		if(!is.null(which.mu))
			mu[which.mu] <- coef(fit)[-(1:ncol(Z))]
		if(verbose) cat("mu", mu, "\n", sep="\t")
		names(fit$df) <- c(uniqueGroups, rep("Offset", length(which.mu)>0))
		if(verbose) cat("df", fit$df,"\n", sep="\t")
		sigma2 = sapply(uniqueGroups, function(i){
					index <- which(groups==i) #& fit$coefficients > beta.thresh
					if(sigma.hat=="BLUP")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ))/(nu + length(index))
					else if(sigma.hat=="df")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ))/(nu + fit$df[i])
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
				sigma2.mu = (nu * sigma0 + sum((mu-0)^2 ) + sum(diag(solve(solve(fit$var)[-(1:ncol(Z)),-(1:ncol(Z))]))))/(nu + length(mu))
			else if(sigma.hat == "REML")
				sigma2.mu = (nu * sigma0 + sum((mu-0)^2 ) + sum(diag(fit$var)[-(1:ncol(Z))]))/(nu + length(mu))
		}
		
		beta = fit$coefficients
		iter = iter+1
	}
	if(iter == max.iter)
		warning("Did not converge after ", max.iter, " iterations.")
	fit$iter[1] <- iter
	fit$sigma2 = sigma0ld
	names(fit$sigma2) <- uniqueGroups
	fit$sigma2.mu = sigma2.mu
	fit$mu = mu
	fit$Z = Z[,order(o)]
	fit$surv = surv
	C <- rbind(diag(1, ncol(Z)),t(as.matrix(MakeInteger(groups)[which.mu]))) ## map from centred to uncentred coefficients 
	fit$groups = groups[order(o)]
	var = fit$var
	var2 = fit$var2
	colnames(var) <- rownames(var) <- colnames(var2) <- rownames(var2) <- rownames(C) <- c(colnames(Z), which.mu)
	colnames(C) <- colnames(Z)
	p <- ncol(Z)
	i <- c(order(o), (1:ncol(var))[-p:-1])
	j <- order(o)
	fit$C <- C[i,j]
	fit$Hinv <- var[i,i] ## Hinv 
	fit$V <- var2[i,i] ## Hinv I Hinv
	fit$z <- (fit$coefficients / sqrt(diag(var)))[i] ## z-scores of centred coefficients
	fit$z2 <- (fit$coefficients / sqrt(diag(var2)))[i] ## z-scores of centred coefficients (var2)
	fit$var = (t(C) %*% var %*% C)[j,j] ## covariance of uncentred coef
	fit$var2 = (t(C) %*% var2 %*% C)[j,j] ## covariance of uncentred coef (var2)
	fit$mu.var = var[-(1:p),-(1:p)] ## covariance of mean
	fit$mu.var2 = var2[-(1:p),-(1:p)] ## covariance of mean (var2)
	fit$means = fit$means[1:p][j]
	fit$coefficients <- (fit$coefficients %*% C)[j]
	names(fit$means) <- names(fit$coefficients) <-  colnames(Z)[j]
	fit$terms <- fit$terms[1:length(uniqueGroups)]
	fit$penalized.loglik <- fit$loglik[2] - fit$penalty[2] - 1/2 * sum(log(fit$sigma2[groups]))
	## Fake call for predict.coxph and survfit.coxph
	call <- match.call()
	if(Z.df){
		call["data"] <- call["Z"]
		formula <- as.formula(paste(as.character(call["surv"]),"~",paste(colnames(Z)[j], collapse="+")))
	}else{
		formula <- as.formula(paste(as.character(call["surv"]),"~",as.character(call["Z"])))
	}
	attr(formula,".Environment") <- parent.frame()
	fit$formula <- formula
	call["formula"] <- call("foo",formula=formula)["formula"]
	fit$terms <- terms(formula)
	fit$call <- call
	class(fit) <- c("CoxRFX", class(fit))
	return(fit)
}

#' Partial risk components
#' @param fit The CoxRFX fit
#' @param newZ New data, defaults to fit$Z
#' @param groups The groups, defaults to fit$groups
#' @return A matrix with the risk per group
#' 
#' @author mg14
#' @export
PartialRisk <- function(fit, newZ=fit$Z, groups=fit$groups) {
	sapply(levels(groups), function(x) {
				ix <- groups == x
				as.matrix(newZ[,ix, drop=FALSE]) %*% coef(fit)[ix] 
			})
}

#' Variance (confidence intervals) of partial risk components
#' @param fit The CoxRFX fit
#' @param newZ New data, defaults to fit$Z
#' @param groups The groups, defaults to fit$groups
#' @param var Variance type. Either 'var' or 'var2'
#' @return A matrix with the confidence interval (prediction variance) for each risk component
#' 
#' @author mg14
#' @export
PartialRiskVar <- function(fit, newZ=fit$Z, groups=fit$groups, var = c("var2","var")) {
	var <- match.arg(var)
	newZ <- newZ - rep(colMeans(newZ), each=nrow(newZ))
	sapply(levels(groups), function(x) {
				ix <- groups == x
				rowSums((as.matrix(newZ[,ix, drop=FALSE]) %*% fit[[var]][ix,ix]) * as.matrix(newZ[,ix, drop=FALSE]))
			})
}

#' Variance components
#' @param fit The CoxRFX fit
#' @param newZ New data, defaults to fit$Z
#' @param groups The groups, defaults to fit$groups
#' @param type Take either the diagonal of the covariance matrix (default), or the rowSums. The latter guaranties that the
#' components sum up to the actual variance, but could be negative in the case of collinearity. The two are equivalent 
#' @param var Which variance estimate to take for the average prediction error. Choices are var2 and var.
#' @return A vector containing the variance components
#' 
#' @author mg14
#' @export
VarianceComponents <- function(fit, newZ = fit$Z[setdiff(1:nrow(fit$Z), fit$na.action),], groups = fit$groups, type = c("diag","rowSums"), var=c("var2","var")){
	var <- match.arg(var)
	risk <- PartialRisk(fit = fit, newZ = newZ, groups = groups)
	type <- match.arg(type)
	#residual <- predict(fit, se.fit=TRUE)$se.fit^2
	newZ <- as.matrix(newZ - rep(colMeans(newZ), each=nrow(newZ)))
	error <- rowSums((newZ %*% fit[[var]]) * newZ)
	
	c <- cov(risk, use="complete")
	if(type=="diag")
		x <- diag(c)
	else
		x <- rowSums(c)
	return(c(x, mean.error=mean(error)))
}

#' Plot variance components
#' @param fit The CoxRFX fit
#' @param col The colors for each component
#' @param groups the groups to be used, if different from the fitted ones.
#' @param type The type of variance compnents: Either 'rowSums' (default) or 'diag'. Rowsums sum up to the actual variance of the linear predictor, but can be
#' negative. Plotting just the 'diag'onal elements guarantees positive components
#' @param conf.int Plot confidence intervals? Default = FALSE
#' @param absolute Whether to plot the absolute variance of the log hazard (TRUE, default), or the relative contribution.
#' @param var Which variance estimated to take for the average prediction error. Choices are var2 and var.
#' @param order Logical or integer. Whether, or if integer how,  the variance components should be ordered.
#' @param digits The number of digits to plot.
  
#' @return NULL
#' 
#' @author mg14
#' @export
PlotVarianceComponents <- function(fit, col=1:nlevels(fit$groups), groups = fit$groups, type="rowSums", conf.int=FALSE, absolute=TRUE, digits=2, var=c("var2","var"), order=TRUE) {
	var <- match.arg(var)
	if(is.null(names(col)))
		names(col) <- levels(groups)
	v <- VarianceComponents(fit, groups=groups, type=type)
	if(is.logical(order))
		if(order==TRUE)
			o <- order(v[levels(groups)], decreasing=TRUE)
		else 
			o <- TRUE
	else 
		o <- order
	v <- v[o]
	if(!absolute)
		v <- v/sum(v)
	vp <- v[v>0]
	vn <- v[v<=0]
	if(conf.int){
		ci <- VarianceComponentsCI(fit, q=c(0.025, 0.975))
		r <- apply(signif(ci[,o], 2), 2, function(x) paste(c("; ", paste(x,collapse="-")), collapse=""))
		rp <- r[v>0]
		rn <- r[v<=0]
	}else{
		ci <- NULL
		rp <- rn <- NULL
	}
	pie(vp, col=col[names(vp)], border=NA, labels=paste(names(vp), " (", round(vp, digits),rp,")",sep=""), radius = sqrt(sum(v)), init.angle=90)
	
	if(length(vn)>0){
		par(new=T)
		pie(c(abs(vn), sum(vp)-sum(vn)), col=c(col[names(vn)],NA), border=NA, labels=paste(names(vn), " (", round(vn, digits),rn,")", sep=""), new=FALSE, density=c(rep(36, length(vn) ),NA), radius = sqrt(sum(v)), init.angle=90)
	}
	invisible(list(vc=v, ci=ci))
}

VarianceComponentsCV <- function(fit, which.coef = grep(":", colnames(fit$Z), invert = TRUE, value=TRUE), type=c("diag","colSums"), method="simple"){
	if(method=="simple"){
		type="diag"
		risk <- PartialRisk(fit = fit)
		Z <- fit$Z
		risk0 <- Z * rep(fit$coef, each=nrow(Z))
		Z <- Z - rep(colMeans(Z), each=nrow(Z))
		#residual <- predict(fit, se.fit=TRUE)$se.fit^2
		V <- fit$var
		colnames(V) <- rownames(V) <- colnames(fit$Z)
		sapply(which.coef, function(i){
					w <- grep(i, colnames(Z), value = TRUE)
					r <- risk
					for(ww in w)
						r[,fit$groups[ww]] <- r[,fit$groups[ww]] - risk0[,ww]
					idx <- !colnames(Z) %in% w
					error <- rowSums((as.matrix(Z[,idx]) %*% V[idx,idx]) * as.matrix(Z[,idx]))
					c <- cov(r, use="complete")
					if(type=="diag")
						x <- diag(c)
					else
						x <- rowSums(c)
					return(c(x, `avg.error`=mean(error)))
				})
	}else{
		sapply(which.coef, function(i){
					w <- grep(i, colnames(fit$Z), value = TRUE, invert=TRUE)
					fit <- CoxRFX(fit$Z[,w], fit$surv, fit$groups[w], sigma0=0.1, nu=0)
					VarianceComponents(fit)
	})}
}


#' Compute concordance for risk components
#' @param fit The CoxRFX fit
#' @param newZ New data, defaults to fit$Z
#' @param groups The groups, defaults to fit$groups
#' @param newSurv The survival object, defaults to fit$surv
#' @return A vector with the concordance of each component
#' 
#' @author mg14
#' @export
PartialC <- function(fit, newZ = fit$Z, newSurv=fit$surv, groups=fit$groups){
	#require(Hmisc)
	c <- sapply(levels(groups), function(x) {
				ix <- groups == x
				r <- as.matrix(newZ[,ix, drop=FALSE]) %*% coef(fit)[ix]
				#rcorr.cens(-r, newSurv)[1]
				c <- survConcordance(newSurv~r)
				c(c$concordance, c$std.err)
			})
	colnames(c) <- levels(groups)
	return(c)
}
	
#' Predict risk with missing data
#' @param fit A CoxRFX fit
#' @param newZ data.frame or matrix with new variables.
#' @param var character. Which variance estimate to use either var = $H^{-1}$, or var2 = $H^{-1}I H^{-1}$. The former is the more conservative choice, the latter (default) seems more accurate, but may underestimate the variance.
#' @param bound logical. Determines whether the imputations should be bound to the range of the variable observed in the original data set.
#' @return A data.frame with columns Expected and Variance 
#' 
#' @author mg14
#' @export
PredictRiskMissing <- function(fit, newZ=fit$Z, var = c("var2","var"), bound=TRUE){
	var <- match.arg(var)
	Sigma <- cov(fit$Z)
	mu <- colMeans(fit$Z)
	rangeZ <- apply(fit$Z,2,range, na.rm=TRUE)
	beta <- fit$coefficients
	newZ <- newZ[,names(beta), drop=FALSE]
	
	.predict <- function(newZ, beta, Sigma, mu, rangeZ){
		missing <- is.na(newZ)
		expectedZ <- newZ
		if(any(missing)){
			if(all(missing)){
				expectedZ[missing] <- mu[missing]
				varianceRisk <- beta[missing] %*% (Sigma[missing,missing]) %*% beta[missing]
			}else{
				s <- Sigma[missing, !missing, drop=FALSE] %*% MASS::ginv(Sigma[!missing, !missing, drop=FALSE])
				expectedZ[missing] <- mu[missing] + s %*% (newZ[!missing] - mu[!missing])
				varianceRisk <- beta[missing] %*% (Sigma[missing,missing] - s %*%  Sigma[!missing, missing] ) %*% beta[missing]
			}
		}else{
			varianceRisk <- 0
		}
		expectedRisk <- expectedZ %*% beta
		if(bound)
			expectedZ <- pmin(pmax(expectedZ, rangeZ[1,]),rangeZ[2,])
		e <- expectedZ - mu
		varianceRisk <- e %*% fit[[var]][1:length(fit$coefficients),1:length(fit$coefficients)] %*% e + varianceRisk
		return(c(expectedRisk, varianceRisk))
	}
	predictions <- t(apply(newZ, 1, .predict, beta, Sigma, mu, rangeZ))
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
#' @param use character Which observations to use for computing the covariance. See cov() for details.
#' @param bound logical. Determines whether the imputations should be bound to the range of the variable observed in the original data set.
#' @return A data.frame of dim(newX) with imputed variables
#' 
#' @author mg14
#' @export
ImputeMissing <- function(X, newX=X, use="pairwise.complete.obs", bound=TRUE){
	Sigma <- cov(X, use=use)
	mu <- colMeans(X, na.rm=TRUE)
	rangeX <- apply(X,2,range, na.rm=TRUE)
	
	l <- ncol(X)
	
	.impute <- function(newX, Sigma, mu, rangeX){
		missing <- is.na(newX)
		expectedX <- newX
		varianceX <- rep(0,l)
		if(any(missing)){
			if(all(missing)){
				expectedX[missing] <- mu[missing]
				varianceX[missing] <- diag(Sigma)
			}else{
			s <- Sigma[missing, !missing, drop=FALSE] %*% MASS::ginv(Sigma[!missing, !missing, drop=FALSE])
			varianceX[missing] <- diag(Sigma[missing, missing] - s %*% Sigma[!missing, missing])
			expectedX[missing] <- mu[missing] + s %*% (newX[!missing] - mu[!missing])
			}
		}
		#return(cbind(expectedX, varianceX))
		if(bound)
			expectedX <- pmin(pmax(expectedX, rangeX[1,]),rangeX[2,])
		expectedX
	}
	imputations <- t(apply(newX, 1, .impute, Sigma, mu, rangeX))
	colnames(imputations) <- colnames(newX) #c("Expected","Variance")
	return(imputations)
}

#' @aliases ImputeMissing
#' @rdname ImputeMissing
#' @export
ImputeXMissing <- function(X, newX=X, use="pairwise.complete.obs", bound=TRUE) .Defunct(ImputeMissing, package = NULL, "ImputeXMissing is now defunct. Please use ImputeMissing instead.")

#' Standardize the magnitude of covariates
#' @param X A matrix or data.frame
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
#' 
#' This plots a CoxRFX model
#' @param x The CoxRFX object
#' @param ... Additional parameters passed to plot(). This can include the non-standard argument order (=1:nlevels(x$groups)), to reorder the groups.
#' @return NULL
#' @method plot CoxRFX
#' @author mg14
#' @export
plot.CoxRFX <- function(x, ...){
	.plt <- function(x, col=c(brewer.pal(9,"Set1"), brewer.pal(8,"Dark2")), order = 1:nlevels(x$groups), xlim=range(coef(x)), xlab="Coefficient",...){
		plot(NA,NA, xlim=xlim, ylim=range(1,1+nlevels(x$groups)), yaxt="n", ylab="",xlab=xlab, ...)
		axis(side=2, at=1:nlevels(x$groups), labels = levels(x$groups)[order], las=2)
		i <- 1
		for(l in levels(x$groups)[order]){
			m <- x$mu[l]
			s <- x$sigma2[l]
			xx <- seq(m-3*sqrt(s),m+3*sqrt(s), l=100)
			y <- i+dnorm(xx, m, sqrt(s))/dnorm(m, m, sqrt(s))*.8
			polygon(xx,y, col=paste(col[order][i],"44", sep=""), border=NA)
			lines(xx, y, col=col[order][i])
			points(x$coefficients[x$groups==l],rep(i, sum(x$groups==l)), col=col[order][i], pch=16, cex=.5)
			lines(rep(m,2),c(0,.8) +i, col=col[order][i])
			i <- i+1
		}
	}
	.plt(x, ...)
}


#' A Wald test for the coefficients of a CoxRFX model
#' 
#' This separately tests the null-hypothesis of being zero on each coefficient in a CoxRFX model using a Wald test. The test
#' statistic is \eqn{z^2 = \beta^2/ Var[\beta]}.
#' 
#' @note Note that there is a lively debate about testing random effects in generalised linear models. See for example
#' http://glmm.wikidot.com/faq
#' @param coxRFX The CoxRFX model
#' @param var Which type of variance estimate to use. The default choice is var2 = H^{-1} I H^{-1}. A more conservative choice is var = H^{-1}.
#' @return A data.frame with columns coef, sd, z and p.value
#' 
#' @author mg14
#' @export
WaldTest <- function(coxRFX, var=c("var2","var")){
	var <- match.arg(var)
	v <- diag(coxRFX[[var]]) 
	z <- coef(coxRFX)/sqrt(v) 
	d <- 1
	p <- pchisq(z^2, d, lower.tail=FALSE)
	data.frame(coef=coef(coxRFX), sd=sqrt(v), z=z, df = d, p.value=p, sig=sig2star(p))
}

#' A summary method for CoxRFX models
#' 
#' This model prints the means and variances for each groups of covariates, as well as the variance components.
#' For the means a Wald test with 1 df is computed testing the null-hypothesis of being zero.
#' 
#' The null-hypothesis of zero variance is tested using a combined Wald test that all coefficients in the group are
#' identical to the mean. Gray (1992) suggests to use \eqn{\beta H^{-1} \beta} as a test statistic in a chi-square
#' test with \eqn{\mathrm{tr}[H^{-1} I]}{tr[H^{-1} I]} df, where H is the Hessian of the penalised model
#' and I is the Hessian of the unpenalised coxph model. Note that all variables taken over the subset of interest only. 
#' As noted by Therneau (2003) this test may be somewhat optimistic. Here, we are using \eqn{z^2 = \sum_i\beta_i^2/H_{ii}},
#' which appears to be a more conservative choice, but the consequences remain to be thoroughly evaluated.  
#' 
#' @note Note that there is a lively debate about testing random effects in generalised linear models. See for example
#' http://glmm.wikidot.com/faq
#' 
#' @references  R. J. Gray (1992). Flexible Methods for Analyzing Survival Data Using Splines, with Applications to Breast Cancer Prognosis. Journal of the American Statistical Association, 87:942-951. http://dx.doi.org/10.1080/01621459.1992.10476248
#' T. M. Therneau, P. M. Grambsch, and V. S. Pankratz (2003). Penalized Survival Models and Frailty. Journal of Computational and Graphical Statistics, 12:156-175. http://dx.doi.org/10.1198/1061860031365
#' @param object A CoxRFX model
#' @param ... Currently unused
#' @return NULL
#' 
#' @author mg14
#' @export
#' @method summary CoxRFX
summary.CoxRFX <- function(object, ...){
	which.mu <- names(object$mu)[object$mu!=0]
	p <- z <- s <- object$mu
	z[which.mu] <- object$mu[which.mu]/sqrt(diag(as.matrix(object$mu.var2)))
	s[which.mu] <- sqrt(diag(as.matrix(object$mu.var2)))
	p <- pchisq(z^2,1,lower.tail=FALSE)
	p[!names(p) %in% which.mu] <- NA
	cat("Means:\n")
	show(format(data.frame(mean=object$mu, sd=s, z=z, p.val=p, sig=sig2star(p)), digits=2))
	cat("\nVariances:\n")
	v <- object$sigma2
	c <- coef(object) - object$mu[object$groups] ## centred coefficients
	chisq <- sapply(split(c^2/diag(object$Hinv)[1:length(c)], object$groups), sum)
	df <- object$df[-(nlevels(object$groups)+1)]
	p <- pchisq(chisq, df, lower.tail=FALSE)
#	f <- as.numeric(table(x$groups)/x$df[-(nlevels(x$groups)+1)])
#	u <- sapply(split(coef(x), x$groups), function(x) sum((x-mean(x))^2)/qchisq(0.025, length(x)))
#	l <- sapply(split(coef(x), x$groups), function(x) sum((x-mean(x))^2)/qchisq(0.975, length(x)))
	show(format(data.frame(sigma2=v, chisq=chisq, df = df, p.val=p, sig=sig2star(p)), digits=2))
	cat("\nPartial log hazard:\n")
	newZ <- object$Z[setdiff(1:nrow(object$Z), object$na.action),]
	p <- PartialRisk(object, newZ = newZ)
	v <- VarianceComponents(object, newZ = newZ)
	e <- colMeans(PartialRiskVar(object, newZ = newZ))
	show(format(data.frame(`Cov[g,g]`=c(diag(cov(p)), TOTAL=NaN), `Sum(Cov[,g])`=c(rowSums(cov(p)),TOTAL=sum(cov(p))), `MSE`=c(e, TOTAL=v[length(v)]), check.names = FALSE),  digits=2))
}

#' Print method for CoxRFX
#' 
#' This function implicitly calls summary.CoxRFX().
#' @param x CoxRFX
#' @param ... Currently unused
#' @return NULL
#' 
#' @author mg14
#' @method print CoxRFX
#' @export
print.CoxRFX <- function(x, ...){
	summary.CoxRFX(x)
}


#' Confidence intervals of variance components
#' 
#' This function numerically computes confidence intervals for the variance components, based sampling on the variances of individual effects.
#' @param fit The CoxRFX model
#' @param newZ The Z (data matrix). Default = fit$Z
#' @param groups The groups. Default fit$groups
#' @param q The quantiles to be evaluated. Default q = c(0.275, 0.5, 0.975), corresponding to the median and symmetric 95\% confidence intervals.
#' @param n The number of samples. Default = 200.
#' @param type The type. Either rowSums or diag(onal).
#' @param absolute Whether the absolute variance of the log hazard should be returned (default), or the relative contributions V_i/sum(V).
#' @param mc.cores The number of mc.cores to use. Default=1.
#' @return A numeric matrix with summary statistics. 
#' 
#' @author mg14
#' @export
VarianceComponentsCI <- function(fit, newZ=fit$Z, groups=fit$groups, q = c(0.025, 0.5, 0.975), type = c("rowSums","diag"), absolute=TRUE, n=200, mc.cores=1){
	type=match.arg(type)
	V <- simplify2array(mclapply(1:n, function(foo){
						set.seed(foo)
						newBeta <- mvtnorm::rmvnorm(1,coef(fit), fit$var2)
						c <- cov(sapply(levels(groups), function(x) {
											ix <- groups == x
											as.matrix(newZ[,ix, drop=FALSE]) %*% newBeta[ix] 
										}))
						if(type=="rowSums")
							rowSums(c)
						else
							diag(c)
					}, mc.cores=mc.cores))
	if(absolute==FALSE)
		V <- V/rep(colSums(V), each=nrow(V))
	apply(V,1, quantile, q)
}