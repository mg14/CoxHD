\name{PartialRiskVar}
\alias{PartialRiskVar}
\title{Variance (confidence intervals) of partial risk components}
\usage{
  PartialRiskVar(fit, newZ = fit$Z, groups = fit$groups,
    var = c("var2", "var"))
}
\arguments{
  \item{fit}{The CoxRFX fit}

  \item{newZ}{New data, defaults to fit$Z}

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

