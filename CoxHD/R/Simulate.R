# Tools for simulating survival in Cox proportional hazards models
# 
# Author: mg14
###############################################################################

#### Basic simulations
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
	FHazInv <- splinefun(c(0,H0$hazard), c(0,H0$time), method="monoH.FC")
	n = length(risk)
	deathTimes = FHazInv(rexp(n, exp(risk))) #predict(hazardDist, rexp(n, exp(risk)))
	
	## Simulate censoring times
	f <- surv
	f[,2] <- 1-f[,2]
	F <- survfit(f~1)
	FCensInv <- splinefun(F$surv, F$time)
	censTimes <- FCensInv(runif(n,0,1)) ## Simulate censoring times
	
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
#' @return A data.frame with extrapolated data
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

