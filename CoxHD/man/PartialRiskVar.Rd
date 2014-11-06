\name{PartialRiskVar}
\alias{PartialRiskVar}
\title{Variance (confidence intervals) of partial risk components}
\usage{
  PartialRiskVar(fit, newX = fit$X, groups = fit$groups,
    var = c("var2", "var"))
}
\arguments{
  \item{fit}{The CoxRFX fit}

  \item{newX}{New data, defaults to fit$X}

  \item{groups}{The groups, defaults to fit$groups}
}
\value{
  A matrix with the confidence interval (prediction
  variance) for each risk component
}
\description{
  Variance (confidence intervals) of partial risk
  components
}
\author{
  mg14
}

