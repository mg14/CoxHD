\name{PartialRisk}
\alias{PartialRisk}
\title{Partial risk components}
\usage{
  PartialRisk(fit, newZ = fit$Z, groups = fit$groups)
}
\arguments{
  \item{fit}{The CoxRFX fit}

  \item{newZ}{New data, defaults to fit$Z}

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

