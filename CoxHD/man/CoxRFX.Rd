\name{CoxRFX}
\alias{CoxRFX}
\title{Cox proportional hazards model with random effects}
\usage{
  CoxRFX(Z, surv, groups = rep(1, ncol(Z)),
    which.mu = unique(groups), tol = 0.001, max.iter = 50,
    sigma0 = 0.1, nu = 0, penalize.mu = FALSE,
    sigma.hat = c("df", "MLE", "REML", "BLUP"),
    verbose = FALSE)
}
\arguments{
  \item{Z}{The data matrix of random effects (n x p)}

  \item{surv}{The survival object (n x 2)}

  \item{groups}{Optional groups as a factor (p) with l
  levels. Default = rep(1, n)}

  \item{which.mu}{Indicator which of the groups should have
  an offset. Default = unique(groups)}

  \item{tol}{The tolerance beyond which to stop}

  \item{max.iter}{The maximal number of iterations}

  \item{sigma0}{The variance of a si-chisq hyperprior on
  the variances.}

  \item{nu}{The df of the variance hyperprior. Default = 0,
  that is no hyperprior.}

  \item{penalize.mu}{Wether to define an N(0,tau)
  hyperprior on the group means.}

  \item{sigma.hat}{Which estimator to use for the
  variances. Default df, other possibilities include MLE,
  REML and BLUP, see details.}

  \item{verbose}{Gives more output.}
}
\value{
  A coxph object with a few extra fields: $groups, $Z,
  $surv, $sigma2 (the variances), $mu (the means), $Hinv
  (the inverse Hessian of the penalised likelihood), $V =
  Hinv I Hinv, the covariance of all coefficients and
  means, $C the map between centred (beta', mu) to beta.
}
\description{
  This function estimates a Cox proportional in which the
  parameters follow normal distributions as discussed by
  Therneau et al. (2003). Multiple groups can be defined
  with different prior mean and variance. The variances of
  the joint distributions are efficiently estimated by an
  EM-type algorithm.
}
\details{
  The values of the means mu_g are estimated using the
  rowSums of Z (within in group) as auxillary variables.

  Different estimators exist for the variances sigma2_g:
  The default is "df", as used by Perperoglou (2014) and
  introduced by Schall (1991). In the M-step of the
  algorithm, this uses sigma^2_g = beta_g beta_g^T/df_g,
  where the degrees of freedom df_g = tr H_{gg} are the
  trace of the Hessian matrix over the elements of group g.
  Alternatives are MLE, REML, and BLUP, as defined by
  Therneau et al. (2003). Simulations indicate that the
  'df' method is most accurate.

  The model is equivalent to coxme(surv ~ (Z1|1) +
  rowSums(Z1) + (Z2|1) + rowSums(Z2) + ...); the coxme
  routine numerically optimises the integrated partial
  likelihood, which may be more accurate, but is
  computationally expensive.
}
\author{
  mg14
}
\references{
  Terry M Therneau, Patricia M Grambsch & V. Shane Pankratz
  (2003) Penalized Survival Models and Frailty, Journal of
  Computational and Graphical Statistics, 12:1, 156-175,
  http://dx.doi.org/10.1198/1061860031365

  A. Perperoglou (2014). Cox models with dynamic ridge
  penalties on time-varying effects of the covariates. Stat
  Med, 33:170-80. http://dx.doi.org/10.1002/sim.5921

  R. Schall (1991). Estimation in generalized linear models
  with random effects. Biometrika, 78:719-727.
  http://dx.doi.org/10.1093/biomet/78.4.719
}

