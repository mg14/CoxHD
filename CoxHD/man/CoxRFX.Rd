\name{CoxRFX}
\alias{CoxRFX}
\title{Cox proportional hazards model with random effects}
\usage{
  CoxRFX(data, surv, groups = rep(1, ncol(data)),
    which.mu = unique(groups), tol = 0.001, max.iter = 50,
    sigma0 = 1, nu = 1, beta.thresh = 0.001,
    penalize.mu = FALSE, sigma.hat = c("df", "p"),
    verbose = FALSE)
}
\arguments{
  \item{data}{}

  \item{surv}{}

  \item{groups}{}

  \item{which.mu}{}

  \item{tol}{}

  \item{max.iter}{}

  \item{sigma0}{}

  \item{nu}{}

  \item{beta.thresh}{}

  \item{penalize.mu}{}

  \item{sigma.hat}{}

  \item{verbose}{}
}
\value{
  A coxph object with a few extra fields.
}
\description{
  Cox proportional hazards model with random effects
}
\author{
  mg14
}

