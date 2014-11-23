\name{VarianceComponents}
\alias{VarianceComponents}
\title{Variance components}
\usage{
  VarianceComponents(fit,
    newZ = fit$Z[setdiff(1:nrow(fit$Z), fit$na.action), ],
    groups = fit$groups, type = c("diag", "rowSums"),
    var = c("var2", "var"))
}
\arguments{
  \item{fit}{The CoxRFX fit}

  \item{newZ}{New data, defaults to fit$Z}

  \item{groups}{The groups, defaults to fit$groups}

  \item{type}{Take either the diagonal of the covariance
  matrix (default), or the rowSums. The latter guaranties
  that the components sum up to the actual variance, but
  could be negative in the case of collinearity. The two
  are equivalent}

  \item{var}{Which variance estimate to take for the
  average prediction error. Choices are var2 and var.}
}
\value{
  A vector containing the variance components
}
\description{
  Variance components
}
\author{
  mg14
}

