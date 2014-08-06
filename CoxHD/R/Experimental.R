# Experimental code, subject to random changes
# 
# Author: mg14
###############################################################################


SimSurvTDNonp <- function(risk, surv, H0 = basehaz(coxph(surv ~ 1)), tplSplit) {
	## Simulate deaths times
	FHazInv <- splinefun(c(0,H0$hazard), c(0,H0$time), method="monoH.FC")
	n = length(risk)
	
	## Simulate censoring times
	nCens <- sum(surv[,3]==0, na.rm=TRUE)
	FSurv <- splinefun(exp(-H0$hazard) ~ H0$time, method="monoH.FC") ## Cum survival dist
	x <- sort(na.omit(surv[surv[,3]==0,1])) ## observed (conditioned) censored times
	t <- table(x) / nCens ## empirical Dist of cond. censoring times
	y <- t / FSurv(unique(x)) ## Need to get unconditional distribution
	FCens <- cumsum(y)
	FCens <- FCens/max(FCens) ## Unconditioned cens.
	FCensInv <- splinefun(unique(x) ~ FCens, method="monoH")
	censTimes <- FCensInv(runif(n,0,1)) ## Simulate censoring times
	
	## TD observations
	tplTimes <- surv[,1]
	tplTimesTT <- FSurv(tplTimes) ## Normalized
	tplIndex <- table(tplSplit) == 2
	
	deathTimes = FHazInv(rexp(n, exp(risk)) + tplTimesTT) ## Predict death conditioned on >= TPL 
	event <- (deathTimes < censTimes)[tlpSplit] + 0
	
	## Put together
	survOut <-  Surv(time = tplTimes, time2 = pmax(0,pmin(deathTimes, censTimes)), event=event)
	return(survOut)
}


