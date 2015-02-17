#' ### Parameters
#' First define some parameters
#+ parameters
set.seed(42)
nParam = 250 # Parameters
nObs = 1000 # Observations
nGroups <- 5
groups <- factor(paste("Group", rep(1:nGroups, each=nParam/nGroups)))

#' Now draw coefficients
#+ coefficients
require(mvtnorm)
mu <- seq(-0.5,0.5,l=nGroups) # Coefficient mean within each group 
sd <- seq(0.1,1, l=nGroups)  # Standard deviations
a <- rnorm(nParam, mean = rep(mu, each=50), sd = rep(sd, each=50)) # Normal coefficients
beta = rbeta(nParam, 1, 20)
Z <- rmvnorm(nObs, mean=rep(0,nParam), sigma = diag(beta) + 1e-3) # Normal covariates with a bit of correlation
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
surv = CoxHD:::SimSurv(risk = risk)
plot(survfit(surv ~1))

#' Maximal concordance
survConcordance(surv ~ risk)

#' ### Fit model
#+ fit
fit <-  CoxRFX(Z, surv, groups = groups) ## takes about 30 s
fit
plot(fit)

#' Coefficients
WaldTest(fit)
plot(a, coef(fit), col=fit$groups)
segments(a, coef(fit) - 2*sqrt(diag(fit$var2)),a, coef(fit) + 2*sqrt(diag(fit$var2)), col=fit$groups)
abline(0,1)

#' Means
plot(mu, fit$mu, col=1:nlevels(groups))
segments(mu,fit$mu - 2*sqrt(diag(fit$mu.var2)) , mu, fit$mu + 2*sqrt(diag(fit$mu.var2)), col=1:nlevels(groups))
abline(0,1)

#' Variances
plot(sd^2, fit$sigma2)
