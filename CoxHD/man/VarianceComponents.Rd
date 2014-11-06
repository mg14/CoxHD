\name{VarianceComponents}
\alias{VarianceComponents}
\title{Variance components}
\usage{
  VarianceComponents(fit, newX = fit$X,
    groups = fit$groups, type = c("diag", "rowSums"),
    var = c("var2", "var"))
}
\arguments{
  \item{fit}{The CoxRFX fit}

  \item{newX}{New data, defaults to fit$X}

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

