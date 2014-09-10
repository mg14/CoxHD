\name{CoxCPSS}
\alias{CoxCPSS}
\title{Complementary pairs stability selection for the Cox proportional hazards model}
\usage{
  CoxCPSS(X, surv, bootstrap.samples = 50, nlambda = 250,
    alpha.weak = 0.5, penalty.factor = rep(1, ncol(X)),
    mc.cores = 1, pi.thr = 0.8,
    control = c("theta", "FDR", p.adjust.methods),
    level = 0.1, simultaneous = TRUE,
    which.error = 1:ncol(X))
}
\arguments{
  \item{X}{matrix or data.frame of the covariats}

  \item{surv}{a Surv() object}

  \item{bootstrap.samples}{Number of bootstrap samples
  (default = 50). Note: These give 100 samples for CPSS.}

  \item{nlambda}{Approximate number of datapoints for the
  LASSO penalty. Default 250.}

  \item{alpha.weak}{Weakness parameter for the stability
  selection. Default 0.5.}

  \item{penalty.factor}{Penalty factor. Allows for reducing
  the penalty on some covariates. Default rep(1,ncol(X)).}

  \item{mc.cores}{Number of cores if parallel computations
  are desired}

  \item{pi.thr}{The stability threshold. Default = 0.8.}

  \item{control}{What type of error control is desired. Use
  "BH" for Benjamini-Hochberg correction.}

  \item{level}{The level for P-value adjustments.}

  \item{simultaneous}{TRUE for complementary pairs
  stability selection. Needed for type-1 error control}

  \item{which.error}{Indeces of covariates used for
  estimating the selection probability under the null.
  Default = 1:ncol(X).}
}
\value{
  CoxCPSS
}
\description{
  Complementary pairs stability selection for the Cox
  proportional hazards model
}
\author{
  mg14
}

