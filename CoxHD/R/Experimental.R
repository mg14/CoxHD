# Experimental code, subject to random changes
# 
# Author: mg14
###############################################################################


SimSurvTDNonp <- function(dataFrame, coef, time0, time1, event, timeTpl, coefTpl) {
	w <- which(timeTpl < time1 & timeTpl > time0)
	index <- c(1:nrow(dataFrame), w) 
	d <- dataFrame[index,]
	d$index <- index
	d$time0 <- c(time0, timeTpl[w])
	d$time1 <- c(pmin(time1, timeTpl, na.rm=TRUE), time1[w])
	d$transplant <- c(rep(0,nrow(dataFrame)), rep(1, length(w)))
	e <- c(event, event[w])
	e[w] <- 0
	d$event <- e
	
	i <- 1
	
	simTime1 <- simTplTimes <- simDeathTimes <- simCensTimes <- rep(NA,nrow(dataFrame))
	#H0 <- basehaz(coxph(Surv(time0, time1, event ) ~ as.matrix(dataFrame) %*% coef), centered = FALSE)
	
	for(idx in list(is.na(timeTpl), !is.na(timeTpl))){
		surv <- Surv(time0[idx], time1[idx], event[idx] )
		
		risk0 <- as.matrix(dataFrame[idx, ]) %*% coef + ifelse(i>1, coefTpl,0) 
		H0 = basehaz(coxph(surv ~ risk0), centered = FALSE)
		
		## Simulate deaths times
		FHazInv <- splinefun(c(0,H0$hazard), c(0,H0$time), method="monoH.FC")
		FHaz <- splinefun(c(0,H0$time), c(0,H0$hazard),  method="monoH.FC")
		n = length(risk0)
		
		## Simulate censoring times
		nCens <- sum(event==0, na.rm=TRUE)
		FSurv <- splinefun(exp(-H0$hazard) ~ H0$time, method="monoH.FC") ## Cum survival dist
		x <- sort(na.omit(time1[event==0])) ## observed (conditioned) censored times
		t <- table(x) / nCens ## empirical Dist of cond. censoring times
		y <- t / FSurv(unique(x)) ## Need to get unconditional distribution
		FCens <- cumsum(y)
		FCens <- FCens/max(FCens) ## Unconditioned cens.
		FCensInvDist <- splinefun(unique(x) ~ FCens, method="monoH")
		simCensTimes[idx] <- FCensInvDist(runif(n,0,1)) ## Simulate censoring times
		
		simDeathTimes[idx] = FHazInv(rexp(n, exp(risk0)))
		#w <- simDeathTimes > simTplTimes & simDeathTimes < simCensTimes
		#w <- sample(seq_along(simDeathTimes), length(tplTimes))
		#simDeathTimes[w] <- FHazInv(FHaz(simDeathTimes[w]) + rexp(length(w), exp(coefTpl + risk0[w])))
		simTime1[idx] <- pmin(simDeathTimes[idx], simCensTimes[idx])

		## TD observations
		if(i > 1){
			w <- which(!is.na(timeTpl))
			tplTimes <- na.omit(timeTpl)
			x <- sort(tplTimes) ## observed (conditioned) censored times
			t <- table(x) / length(tplTimes) ## empirical Dist of cond. censoring times
			#y <- t / FSurv(unique(x)) ## Need to get unconditional distribution
			FTpl <- cumsum(t)
			FTpl <- FTpl/max(FTpl) ## Unconditioned cens.
			FTplDist <- splinefun(FTpl ~ unique(x), method="monoH")
			FTplInvDist <- splinefun(unique(x) ~ FTpl, method="monoH")
			#simTplTimes <- FTplInvDist(runif(n,0,1)) ## Simulate censoring times
			
			#simTplTimes[w] <- pmin(FTplInvDist(runif(length(w),0,pmax(0,pmin(1,FTplDist(simTime1[w]))))), simTime1[w]-0.5)
			simTplTimes[w] <- runif(length(w), min(simTime1[w]), simTime1[w])
		}
		i <- i+1
	}
	## Put together
	index <- c(1:nrow(dataFrame), w) 
	d <- dataFrame[index,, drop=FALSE]
	d$index <- index
	d$time0 <- c(time0, simTplTimes[w])
	d$time1 <- c(pmin(simTime1, simTplTimes, na.rm=TRUE), simTime1[w])
	d$transplant <- c(rep(0,nrow(dataFrame)), rep(1, length(w)))
	evt <- simDeathTimes < simCensTimes
	e <- c(evt, evt[w])
	e[w] <- 0
	d$event <- e
	return(d)
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

GetPairs <- function(names, scope){
	pairs <- strsplit(scope, "\\.")
	lapply(pairs, function(s){
				a <- grep(paste(s[1],"$", sep=""), names, value=TRUE)
				b <- grep(paste(s[2],"$", sep=""), names, value=TRUE)
				c(a,b)
			})
}


WaldTest <- function(coxRFX, var=c("var2","var")){
	var <- match.arg(var)
	v <- diag(coxRFX[[var]]) 
	z <- coef(coxRFX)/sqrt(v) 
	d <- 1
	p <- pchisq(z^2, d, lower.tail=FALSE)
	data.frame(coef=coef(coxRFX), sd=sqrt(v), z=z, df = d, p=p, sig=sig2star(p))
}

show.CoxRFX <- function(x){
	which.mu <- names(x$mu)[x$mu!=0]
	p <- z <- s <- x$mu
	z[which.mu] <- x$mu[which.mu]/sqrt(diag(x$mu.var2))
	s[which.mu] <- sqrt(diag(x$mu.var2))
	p <- pchisq(z^2,1,lower.tail=FALSE)
	p[!names(p) %in% which.mu] <- NA
	cat("Means:\n")
	show(format(data.frame(mean=x$mu, sd=s, z=z, p.val=p, sig=sig2star(p)), digits=2))
	cat("\nVariances:\n")
	v <- x$sigma2
	c <- coef(x) - x$mu[x$groups] ## centred coefficients
	chisq <- sapply(split(c^2/diag(x$Hinv)[1:length(c)], x$groups), sum)
	df <- x$df[-(nlevels(x$groups)+1)]
	p <- pchisq(chisq, df, lower.tail=FALSE)
#	f <- as.numeric(table(x$groups)/x$df[-(nlevels(x$groups)+1)])
#	u <- sapply(split(coef(x), x$groups), function(x) sum((x-mean(x))^2)/qchisq(0.025, length(x)))
#	l <- sapply(split(coef(x), x$groups), function(x) sum((x-mean(x))^2)/qchisq(0.975, length(x)))
	show(format(data.frame(sigma2=v, chisq=chisq, df = df, p.val=p, sig=sig2star(p)), digits=2))
}

print.CoxRFX <- function(x){
	show.CoxRFX(x)
}

concordanceFromVariance <- function(x) {
	.cfv <- function(x){
		stopifnot(x >= 0)
		if(x == 0)
			return(.5)			
		f <- function(x, sigma2) dnorm(x, sd=sqrt(2*sigma2)) / (1 +  exp(- x))
		2*integrate(f, 0, 10*sqrt(x), sigma2=x)$value
	}
	sapply(x, .cfv)
}

SimCoef <- function(coxRFX=NULL, groups = coxRFX$groups, mu=coxRFX$mu, sigma2=coxRFX$sigma2){
	c <- rnorm(length(groups), mu[groups], sqrt(sigma2[groups]))
	names(c) <- names(groups)
	c
}