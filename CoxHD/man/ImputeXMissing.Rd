\name{ImputeXMissing}
\alias{ImputeXMissing}
\title{Covariance-based imputation of missing variables}
\usage{
  ImputeXMissing(X, newX)
}
\arguments{
  \item{X}{orignial data set}

  \item{newX}{The data.frame of covariates}
}
\value{
  A data.frame of dim(newX) with imputed variables
}
\description{
  This functinon imputes missing variables based on the
  covariance of a previous data set. This can be useful if
  a CoxRFX model (or any other) has been fit on a large
  data set with correlated variables and one wants to use
  this model on a different data sets with fewer or
  randomly missing covariates.
}
\author{
  mg14
}

