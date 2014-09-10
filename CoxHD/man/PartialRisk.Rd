\name{PartialRisk}
\alias{PartialRisk}
\title{Partial risk components}
\usage{
  PartialRisk(fit, newX = fit$X, groups = fit$groups)
}
\arguments{
  \item{fit}{The CoxRFX fit}

  \item{newX}{New data, defaults to fit$X}

  \item{groups}{The groups, defaults to fit$groups}
}
\value{
  A matrix with the risk per group
}
\description{
  Partial risk components
}
\author{
  mg14
}

