#' The Cox RFX model
#' =========================

#+ echo=FALSE, cache=FALSE
opts_chunk$set(cache=TRUE, autodep=TRUE)
options(width=120)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
		})
opts_chunk$set(dev=c('png','pdf'), fig.ext=c('png','pdf'), fig.width=4, fig.height=4, smallMar=TRUE)


#' ### Preliminaries
#+ preliminaries, cache=FALSE
set.seed(42)
library(CoxHD)
library(Hmisc)
library(RColorBrewer)
library(mvtnorm)

#' ## A random effects model
#' To estimate survival we use an extension of the Cox proportional hazards model. In this model the hazard, given a set of covariates $Z$
#' is expressed as
#' $$\lambda = \lambda_0(t) \exp(Z u) $$
#' where $u$ denotes the risk coefficients.The 
#' random effects model assumes that the coefficients $u_j$ are from a **shared distribution**
#' $$u_j \sim N(\mu, \sigma^2).$$
#' This assumption make the model a hierarchical model, which is often referred to as a random effects model. The assumption of a shared distribution with a finite variance effectively regularises
#' the model and allows for estimating parameters even in  high-dimensional cases. 
#'  
#' A second idea of our model is that our the set of covariates $Z = (Z_1,...,Z_g)$ can be split into $g$ **groups**, each with a different distribution of effect sizes, 
#' $$ u_j \sim N(\mu_g(j), \sigma^2_g(j))$$.
#' 
#' The means $\mu_g$ can be estimated utilising the decomposition of the log-hazard, $h_i = \sum_g h_{ig}$, $h_{ig} = \sum_{j\in g} X_{ij}u_j =  \sum_{j\in g} X_{ij}(u_j - \mu_g) + \sum_{j\in g} X_{ij}\mu$. 
#' So one can estimate each $\mu_g$ as the coefficient of the introducing auxillary variables $\bar X_g = \sum_{j\in g} X_{ij}$, which are simply the row-wise sums of the covariates of each group. 
#' 
#' To estimate the centred parameters $u'_j = u_j - \mu_{g(j)}$, we exploit that the constraint induced by the normal prior distributions is equivalent to a ridge penalty. 
#' This gives the MAP estimates
#' $$u_j^* \mid \sigma_g(j)^2, Z = u_j^{ridge},$$
#' which are implemented in R using the ridge() function in coxph().
#' 
#' To estimate the variances $\sigma_g^2$ we can then iterate between the MAP estimates $u^*$ and
#' $$\sigma_g^2 \mid u_j^*,j\in g = \frac{\sum_ju_j^2 }{df} $$
#' where $df$ are the effective degrees of freedom, $df = \mathrm{tr} [(H_{gg}+\sigma_g^2 I)(H_{gg})^{-1}]$, $H$ being the Hessian of the total likelihood, evaluated for elements of group $g$.
#' 
#' The covariance of the coefficients $u_j = u_j' + \mu_{g(j)}$ can be derived from the $(p+|g|) \times (p+|g|)$ covariance matrix $V$ of the coeffient vector $(u', \mu)$, 
#' $\Sigma_u = B^T V B$, where $B = (\mathbb{1}_p, I(j\in g)$ maps $(u', \mu)$ to $u$. There exist two different estimates of $V$, either $V=H^{-1}$ or  $V = H^{-1} \mathcal{I} H^{-1}$, but 
#' we find that the former tends to be overly conservative.
#' 
#' Optionally one may define a hyperprior for $\sigma^2_g \sim \operatorname{si}\chi^2(\nu, \sigma_0^2)$, which can help stabilize the estimates.
#' 
#' The implementation of the model is straightforward: 
CoxRFX

#' ## Simulations
#' First define some parameters
#+ parameters
nParam = 250 # Parameters
nObs = 1000 # Observations
nGroups <- 5
groups <- factor(paste("Group", rep(1:nGroups, each=nParam/nGroups)))

#' Now draw coefficients
#+ coefficients, cache=FALSE
mu <- seq(-0.5,0.5,l=nGroups) # Coefficient mean within each group 
sd <- seq(0.1,1, l=nGroups)  # Standard deviations
a <- rnorm(nParam, mean = rep(mu, each=50), sd = rep(sd, each=50)) # Normal coefficients
beta = rbeta(nParam, 1, 20)
#Z = sapply(beta, function(x) rnorm(n=nObs, mean=0, sd = sqrt(x))) # Normal covariates
Z <- rmvnorm(nObs, mean=rep(0,nParam), sigma = diag(beta) + 1e-3)
Z[] <- Z > quantile(Z, 0.75) ## Make binary
cor(Z[,1:5])

#' Simulated risk
#+ risk
risk = Z %*% a
a <- a / sd(risk) # standardize
mu <- mu / sd(risk) # standardize
sd <- sd / sd(risk) # standardize
risk <- risk / sd(risk) # standardize
head(risk)

#' By group
#+ groups
riskComponents <- sapply(levels(groups), function(g)  Z[,groups==g] %*% a[groups==g])
cov(riskComponents)
rowSums(cov(riskComponents))

#' ### Simulate survival
#+ survival
CoxHD:::SimSurv
surv = CoxHD:::SimSurv(risk = risk)
plot(survfit(surv ~1))

#' Maximal concordance
survConcordance(surv ~ risk)

#' ### Fit model
#+ fit
fit = CoxRFX(Z, surv, groups = groups)

fits <- lapply(1:10, function(x){
			surv = CoxHD:::SimSurv(risk = risk)
			CoxRFX(Z, surv, groups = groups, sigma0 = 0.1, nu=0)
		})


#' Plot estimates
#+ plotEstimates, fig.height=8, fig.width=8
par(mfrow=c(2,2))
boxplot(coef(fit) ~ groups)

plot(rep(sd,10), sqrt(sapply(fits, `[[`, "sigma2")), xlab='sd', ylab="estimate")
abline(0,1)

plot(rep(mu,10), sapply(fits, `[[`, "mu"), xlab='mu', ylab="estimate")
abline(0,1)

plot(a, coef(fit),  col= brewer.pal(5,"Set1")[groups])
abline(0,1)

#' Compared to standard coxph
#+ fig.height=4, fig.width=4
coxfit <- coxph(surv ~ Z)
plot(a, coef(coxfit), pch="")
abline(0,1)
arrows(a,coef(coxfit),a, coef(fit), length=0.1)
mean((a-coef(fit))^2)
mean((a-coef(coxfit))^2)

#' Risk contributions
#+ riskContributions, fig.height=4, fig.width=4
estRiskComponents <- PartialRisk(fit)
estRisk <- rowSums(estRiskComponents)
plot(estRisk, risk)
plot(survfit(surv ~ cut(estRisk, quantile(estRisk, seq(0,1,l=4)))))

#' ### Variance components
#' The variance of the log hazard $h = Z u$ can be written as
#' \[
#' \begin{aligned}
#' Var[h] &= E[Var[h|u]] + Var[E[h|u]] \cr
#' &= E_u Var_Z [Zu] +  Var_u[E[Z]u] \cr
#' &\approx \hat u^T \Sigma_X \hat u + E[Z]^T \Sigma_u E[Z]\cr 
#' &= V + E[\sigma^2_0].
#' \end{aligned}
#' \]
#' The interpretation of the two terms are the following: The first is the predicted variance $V$ using fixed parameters $\hat u$ stemming from the variance in the covariates. 
#' The second term is the uncertainty of the prediction resulting from the variance in the coefficient estimates. 
#' It is tempting to approximate the moments of $u$ with the estimates from the fit. The covariance $\Sigma_Z$ can be estimated across the rows of $Z$. 
#' The covariance of $u$ can be estimated in the fitting process.
#' \[
#' \begin{aligned}
#' Var[h] & \approx Var_Z [Z^T \hat u] + E[Z \hat\Sigma_u Z] \cr
#' &= \hat u^T \Sigma_Z \hat u + E[Z^T \hat\Sigma_u Z] 
#' \end{aligned}
#' \]
#' The explained variance $V$ involves summations over all terms $j$ of $Z_ij$. These can be partitioned into groups g,
#' \[
#' \begin{aligned}
#' V &= \sum_j \sum_k u_j u_k {\Sigma_Z}_{jk} \cr
#' &= \sum_g \sum_h V{gh}, 
#' \end{aligned}
#' \]
#' where $V{gh} = \sum_{j \in g}\sum_{k \in h} u_j u_k ({\Sigma_Z})_{jk}$. $V_{gh}$ is the covariance matrix of the risk contributions. 
#' Ignoring off-diagonal terms in $V_{gh}$ the risk can hence be written as
#' $$ Var[h] \approx \sum_g V_{gg} + E[\sigma^2_0].$$
#' The actual implementation is straightforward:
VarianceComponents
varComp <- VarianceComponents(fit)
pie(varComp,main = "Variance components")

#+ dotchart, fig.height=5, fig.width=5
comp <- lapply(1:10, function(i){
			surv <-  CoxHD:::SimSurv(risk = risk)
			fit <-  CoxRFX(Z, surv, groups = groups, sigma0 = 0.1, nu=0)
			riskComponents <- sapply(levels(groups), function(g)  Z[,groups==g] %*% a[groups==g])
			varComp <- VarianceComponents(fit)
			c(varComp, Residual=mean((predict(fit) - risk)^2))
		})

dotchart(c(diag(cov(riskComponents)), Error=mean((predict(fit) - risk)^2)), pch=19, main = "Variance components", xlim=c(0,max(varComp)),)
for(v in comp)
	points(v[1:6],seq_along(varComp))
legend("topright", pch=c(19,1), c("True","Estd"), bty="n")

#' ### Comparison to frailty models
#' For a factorial set of covariates (i.e., $Z_{ij} \in \{0,1\}; \sum_j Z_{i,j}=1 \forall i$), the model is equivalent to a frailty model
#+ frailty, fig.height=8, fig.width=8
set.seed(42)
par(mfrow=c(3,3))
for(nLevels in c(5,10,50))
	for(n in c(2,10,50)){
		f <- factor(rep(1:nLevels, each=n))
		Z1 <- sapply(1:nLevels, `==`, f) +0
		b <- rnorm(nLevels)
		r <- Z1 %*% b
		s <- CoxHD:::SimSurv(r)
		frailtyFit <- coxph(s ~ frailty(f, "gaussian"))
		rfxFit <- CoxRFX(Z1, s, sigma0=1, which.mu=NULL)
		if(nLevels<=5)
			y <- frailtyFit$coef
		else
			y <- frailtyFit$frail
		plot(coef(rfxFit), y, xlab='CoxRFX est', ylab="Frailty est", pch=19, xlim=range(coef(rfxFit)), ylim=range(coef(rfxFit)))
		title(main=paste(nLevels, "levels,", n,"obs"))
		points(coef(rfxFit), b, pch=1)
		abline(0,1)
	}

#' ### Comparison to coxme
#+ coxme
library(coxme)
ZZ <- lapply(levels(groups), function(l) as.matrix(Z[,groups==l]))
system.time(fitme <- coxme(surv ~ (ZZ[[1]]|1) + (ZZ[[2]]|1) + (ZZ[[3]]|1) +(ZZ[[4]]|1) +(ZZ[[5]]|1) + sapply(ZZ, rowSums)))
system.time(fit <- CoxRFX(Z, surv, groups))

#' CoxRFX is about 100x faster.
#' The estimates of the coefficients are almost identical:
plot(fit$coefficients - fit$mu[fit$groups], unlist(fitme$frail), xlab="coef. CoxRFX",ylab="coef. coxme")
abline(0,1)

#' The variance estimate of CoxRFX seem to be slightly larger:
plot(fit$sigma2, unlist(fitme$vcoef), xlab="sigma2 CoxRFX",ylab="sigma2 coxme")
abline(0,1)


#' ## Missing data
#' If we a proportion missing data $X_m$ and observed data $X_o$, we can impute predictions using the correlation structure $\Sigma_X$ of the training data set.
#' Under normality assumptions one can use the covariance matrix from the training data to obtain the conditional (posterior) distribution
#' of the missing data:
#' \[ \begin{aligned} 
#' X_m | X_o & \sim N(\mu_m^*, \Sigma_m^*) \cr
#' \mu_m^* &= \mu_m + \Sigma_{mo} \Sigma_{oo}^{-1} (X_o - \mu_o) \cr
#' \Sigma_m^* &= \Sigma_{mm} - \Sigma_{mo} \Sigma_{oo}^{-1} \Sigma_{om}.
#' \end{aligned} \]
#' 
#' As the risk predictions are linear in $X$ the distribution of the prediction error is given by
#' \[ \begin{aligned}
#' E[\lambda | X_o, u] &= X_o u_o + \mu_m^* u_m \cr
#' Var[\lambda | X_o, u ] &= u_m^T \Sigma_m^* u_m
#' \end{aligned} \]
#' The uncertainty in $u$ introduces a second term 
#' \[ Var[\lambda | X_o] = Var[\lambda | X_o, u ] + (X_o - \mu_o, \mu_m^* - \mu_m)^T \Sigma_u (X_o -\mu_o, \mu^*_m -\mu_m). \]
#' The implementation is
#+ predictRiskMissing
PredictRiskMissing

#' If nothing is missing this is equivalent to the standard variance estimates of the glm predictors $\sigma^2_i = X_ij {\Sigma_u}_{jk} X_ik$:
#+ predict, fig.height=4, fig.width=4
plot(PredictRiskMissing(fit, var="var")[,2], predict(fit, se.fit=TRUE)$se.fit^2, xlab = "sigma^2 CoxHD", ylab="sigma^2 coxph")

#+ missing, fig.height=4, fig.width=8, fig.show="asis"
for(m in c(0, 0.05, 0.1, 0.5, 0.75, 0.9, 0.95)){
	par(bty="L", mfrow=c(1,2))
	nNewZ <- 500
	newZ <- fit$Z[1:nNewZ,]
	newZ[sample(1:length(newZ),length(newZ) * m, replace = FALSE)] <- NA
	p <- PredictRiskMissing(fit, newZ[1:nNewZ,])
	p0 <- risk[1:nNewZ]
	plot(p0[1:100], p[1:100,1], xlab="True risk", ylab= "Predicted risk", main=paste(100*m,"% missing", sep=""), ylim=range(p0), xlim=range(p0))
	segments(p0[1:100], p[1:100,1] - 2*sqrt(p[1:100,2]),p0[1:100], p[1:100,1] + 2*sqrt(p[1:100,2]), col="#88888844", lwd=4 )
	abline(0,1)
	legend("topleft", bty="n", c(paste("Cor =", round(cor(p[,1],p0),2), collapse=""), paste("RSS =", round(mean((p0 - p[,1])^2), 2), collapse=""), paste("^RSS =", round(mean(p[,2]), 2), collapse="") ))
	plot(survfit(surv[1:nNewZ] ~ (p[,1] >  median(p[,1]))), xlab="Time", ylab="Survival")
	legend("topright", bty="n",paste("C = ", round(rcorr.cens(-p[,1], surv[1:nNewZ])[1],2), collapse=""))
}

#+ rho, fig.height=4, fig.width=10, fig.show="asis"
library(Matrix)
for(rho in -4:-1){
	m <- 0.75
	par(bty="L", mfrow=c(1,3))
	Z <- rmvnorm(nObs, mean=rep(0,nParam), sigma = diag(beta) + as.matrix(bdiag(lapply(1:nGroups, function(i) matrix(10^rho, nrow=nParam/nGroups, ncol=nParam/nGroups)))))
	Z[] <- Z > quantile(Z, 0.75)
	risk = Z %*% a
	risk <- risk / sd(risk)
	s <- CoxHD:::SimSurv(risk)
	fit <- CoxRFX(Z, s, groups=groups)
	nNewZ <- 500
	newZ <- fit$Z[1:nNewZ,]
	newZ[sample(1:length(newZ),length(newZ) * m, replace = FALSE)] <- NA
	p <- PredictRiskMissing(fit, newZ[1:nNewZ,], var="var2")
	p0 <- risk[1:nNewZ]
	plot(p0[1:100], p[1:100,1], xlab="True risk", ylab= "Predicted risk", main=paste(100*m,"% missing, rho=", 10^rho, sep=""), ylim=range(p0), xlim=range(p0))
	segments(p0[1:100], p[1:100,1] - 2*sqrt(p[1:100,2]),p0[1:100], p[1:100,1] + 2*sqrt(p[1:100,2]), col="#88888844", lwd=4 )
	abline(0,1)
	legend("topleft", bty="n", c(paste("Cor =", round(cor(p[,1],p0),2), collapse=""), paste("RSS =", round(mean((p0 - p[,1])^2), 2), collapse=""), paste("^RSS =", round(mean(p[,2]), 2), collapse="") ))
	plot(a/sd( Z %*% a), fit$coefficients, xlab="True coef", ylab="Est. coef")
	abline(0,1)
	plot(survfit(s[1:nNewZ] ~ (p[,1] >  median(p[,1]))), xlab="Time", ylab="Survival")
	legend("topright", bty="n",paste("C = ", round(rcorr.cens(-p[,1], s[1:nNewZ])[1],2), collapse=""))
}

#' # TODO
#' * Mixed model implementation.