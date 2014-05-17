# Complementary pairs stability selection for the Cox model
# 
# Author: mg14
###############################################################################




#### Stability selection
#' Complementary pairs stability selection for the Cox proportional hazards model
#' @param X matrix or data.frame of the covariats
#' @param surv a Surv() object
#' @param bootstrap.samples Number of bootstrap samples (default = 50). Note: These give 100 samples for CPSS.
#' @param nlambda Approximate number of datapoints for the LASSO penalty. Default 250.
#' @param alpha.weak Weakness parameter for the stability selection. Default 0.5.
#' @param penalty.factor Penalty factor. Allows for reducing the penalty on some covariates. Default rep(1,ncol(X)).
#' @param mc.cores Number of cores if parallel computations are desired
#' @param pi.thr The stability threshold. Default = 0.8. 
#' @param control What type of error control is desired. Use "BH" for Benjamini-Hochberg correction.
#' @param level The level for P-value adjustments.
#' @param simultaneous TRUE for complementary pairs stability selection. Needed for type-1 error control
#' @param which.error Indeces of covariates used for estimating the selection probability under the null. Default = 1:ncol(X).
#' @return CoxCPSS
#' 
#' @author mg14
#' @export
CoxCPSS <- function(X, surv, bootstrap.samples=50, nlambda=250, alpha.weak=0.5, penalty.factor = rep(1,ncol(X)), mc.cores=1, pi.thr=0.8, control = c("theta","FDR", p.adjust.methods), level=0.1,  simultaneous = TRUE, which.error = 1:ncol(X)) {
	d = floor(nrow(X)/2)
	control <- match.arg(control)
	
	#### GLMNET dry run
	alpha.net <- 1 ## 1: LASSO, <1 : elastic net
	
	path <- glmnet(x=as.matrix(X[1:d,]), y = as.matrix(data.frame(time=surv[1:d,1] , status = surv[1:d,2])), family="cox", standardize = F, alpha = alpha.net, nlambda = nlambda, penalty.factor = penalty.factor)	
	
	lambda <- path$lambda    #10^(seq(-3,0, length.out=500))
	nlambda <- length(lambda)
	
	gc()
	x <- as.matrix(X) - rep(colMeans(X), each=nrow(X))
	
	l <- mclapply(1:bootstrap.samples, function(b){
				s <- sample(nrow(x), d, replace=F)
				i <- sample(3,1)
				w <- runif(ncol(x),alpha.weak,1)
				#x / rep(pmax(apply(x,2,max),1), each=nrow(x))
				x <- jitter(x) * rep(w, each=nrow(x))
				bath <- try(glmnet(x=x[s,], y = as.matrix(data.frame(time = surv[s,1], status = surv[s,2])), lambda=lambda, family="cox", standardize=F, alpha = alpha.net, nlambda=250,  penalty.factor = penalty.factor))
				res = list()
				if(class(bath)[1] != "try-error" & all(dim(coef(bath))==c(ncol(x), nlambda))){
					res[[1]] = coef(bath) != 0
					rownames(res[[1]]) = colnames(X)
					colnames(res[[1]]) = NULL
				}
				else{
					cat("x")					
					return(NULL)
				}
				if(simultaneous == TRUE){
					bath <- try(glmnet(x=x[-s,], y = as.matrix(data.frame(time = surv[-s,1], status = surv[-s,2])), lambda=lambda, family="cox", standardize=F, alpha = alpha.net, nlambda=250,  penalty.factor = penalty.factor))
					if(class(bath)[1] != "try-error" & all(dim(coef(bath))==c(ncol(x), nlambda)))
					{
						res[[2]] = coef(bath, s=lambda) != 0
						rownames(res[[2]]) = colnames(X)
						colnames(res[[2]]) = NULL
					}
					else{
						cat("X")
						return(NULL)
					}
				}
				cat(".")
				return(res)
			}, mc.cores = mc.cores)
	cat("\n")
	Pr <- matrix(0, dim(path$beta)[1], length(lambda))
	if(simultaneous) Pr_sim <- Pr
	else Pr_sim = NULL
	M <- 0
	for(b in l)
		if(!is.null(b)){
			for(bb in b) 
			{ 
				Pr <- Pr + as.matrix(bb) 
				M <- M+1
			}
			if(length(b) == 2)
			{
				Pr_sim = Pr_sim + as.matrix (b[[1]] & b[[2]])
			}
		}
	B = length(na.omit(l))
	if(M!= (1+simultaneous) * bootstrap.samples) warning("Number of samples smaller than number of bootstrap samples. Some have failed")
	Pr <- Pr/M	
	Pr_sim = Pr_sim/B
	Lambda <- NULL
	res = list(X=X, surv=surv, lambda=lambda, Lambda=Lambda, Pr=Pr, Pi=apply(Pr, 1, max), bootstrap.samples=B, n=nrow(X), pi.thr=pi.thr, level=level, control=control, Pr_sim=Pr_sim, M=M,  penalty.factor =  penalty.factor , alpha.weak = alpha.weak, P=NULL, which.error=which.error)
	class(res) = "CoxCPSS"

	ErrorControlCPSS(res, control = control, level = level)
}

#' Type-1 error control for CoxCPSS
#' 
#' This function implements P-value based error control for the CPSS (Shah and Samworth, 2013). This transforms the selection probability into P-values
#' for L-concave distributions. Common multiple testing adjustments can then be performed on those P-values
#' @param coxCPSS A CoxCPSS object
#' @param control Which type of control
#' @param level The level for P-value adjustments
#' @param pi.thr The threshold for stable variables
#' @param which.error Which covariates are used for estimating the parameter of the null-distribution
#' @return CoxCPSS
#' 
#' @author mg14
#' @export
ErrorControlCPSS <- function(coxCPSS,  control = coxCPSS$control, level=coxCPSS$level, pi.thr=coxCPSS$pi.thr, which.error=coxCPSS$which.error) {
	Pr <- coxCPSS$Pr
	M <- coxCPSS$M
	
	q = colSums(Pr)
	errorControlIdx <- coxCPSS$penalty.factor>0 & 1:ncol(coxCPSS$X) %in% which.error
	theta <- colMeans(Pr[errorControlIdx, ])
	#theta <- apply(Pr[stabCox$penalty.factor>0, ],2,median)
	
	
	if(control == "FDR"){
		if(pi.thr <= 3/4)
			coxCPSS$Lambda = which(q/nrow(Pr) * 1 / (2 * (2*pi.thr - 1 - 1/(2*M) )) < level ) #TODO: double check
		else
			coxCPSS$Lambda = which(q/nrow(Pr) * 4 * (1 - pi.thr + 1/(2*M)) / (1 + 1/M ) < level )
	}else if(control == "theta"){
		coxCPSS$Lambda <- which(theta <= level)
	}else if(control %in% p.adjust.methods){
		P <- sapply(seq(0,round(max(theta),2),.01), function(t){
					if(t > 0 & t < .49) ## TODO: why not .5?
						c(1,minD(t, M/2))
					else if( t >= .49)
						rep(1,M+1)
					else
						c(1,rep(0, M))
				})
		coxCPSS$Lambda <- which(sapply(1:ncol(Pr), function(i){
							pr <- Pr[,i]
							theta <- round(theta[i],2)
							all(which(pr[coxCPSS$penalty.factor>0] > pi.thr) %in% which(p.adjust(P[pr[coxCPSS$penalty.factor>0] * M +1,theta*100 +1], method=control) < level))
						}))
		coxCPSS$P <- P
		i <- max(coxCPSS$Lambda)
		pr <- Pr[,i]
		t <- round(theta[i],2)
		coxCPSS$Pval <- P[pr * M +1,t*100 +1]
		names(coxCPSS$Pval) <- rownames(Pr)
	}
	coxCPSS$Pi <- apply(Pr[,coxCPSS$Lambda], 1, max)
	coxCPSS$control <- control
	coxCPSS$level <- level
	coxCPSS$pi.thr <- pi.thr
	coxCPSS$theta <- theta
	coxCPSS
}

#' Plot a CoxCPSS model
#' @param x a CoxCPSS model 
#' @param xlab 
#' @param ylab 
#' @param lty 
#' @param log 
#' @param col 
#' @param xlim 
#' @param ylim 
#' @param ... 
#' @return NULL
#' @S3method plot CoxCPSS
#' 
#' @author mg14
#' @export
plot.CoxCPSS = function(x, xlab='1/lambda', ylab="Selection probability", lty = rep(1, ncol(x$X)), log="", col = c("black","grey"), xlim=range(1/x$lambda[x$Lambda]), ylim=c(0,1),...){
	#attach(stacoxph)
	plot(NA,NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, log = log, ... )
	for(i in 1:nrow(x$Pr)){
		lines(1/x$lambda, x$Pr[i,], lwd = 1, col=ifelse(x$Pi[i] > x$pi.thr, col[1], col[2]), lty = lty[i], ...)
	}
	lines(1/x$lambda, x$theta, lwd = 1, lty= 3)
	#mtext(at=Pr[names(Pi),][Pi > pi_thr, max(which(Lambda)],.01), names(Pi)[Pi> pi_thr],side=4, cex=.75)
}

#' CoxCPSS with interaction terms
#' @param X 
#' @param surv 
#' @param scope 
#' @param ... 
#' @return CoxCPSS
#' 
#' @author mg14
#' @export
CoxCPSSInteractions <- function(X, surv, scope = 1:ncol(X),...){
	#R <- apply(X, 2, sample)
	#s1 <- StabCox(cbind(X,R=R), surv, control="BH", which.error = ncol(X) + 1:ncol(X), ...)
	s1 <- CoxCPSS(X, surv, control="BH")
	w <- which(s1$Pi > s1$pi.thr)
	i <- intersect(scope, w)
	I <- MakeInteractions(X[,i],X[,i])[,as.vector(upper.tri(matrix(0,ncol=length(i), nrow=length(i))))]
	I <- I[,colSums(I) > 0]
	Z <- cbind(X, I)
	#R <- apply(Z, 2, sample)
	#Z <- cbind(Z, R=R)
	penalty <- rep(1, ncol(Z))
	penalty[i] <- 0
	#StabCox(Z, surv, penalty.factor = penalty, control = "BH", which.error = (ncol(X) + ncol(I) + 1):ncol(Z) , ...)
	CoxCPSS(Z, surv, penalty.factor = penalty, control = "BH", ...)
}

