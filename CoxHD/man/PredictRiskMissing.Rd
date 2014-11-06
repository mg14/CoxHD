\name{PredictRiskMissing}
\alias{PredictRiskMissing}
\title{Predict risk with missing data}
\usage{
  PredictRiskMissing(fit, newX = fit$X,
    var = c("var2", "var"))
}
\arguments{
  \item{fit}{A CoxRFX fit}

  \item{newX}{}

  \item{var}{Which variance estimate to use either var =
  $H^{-1}$, or var2 = $H^{-1}I H^{-1}$. The former is the
  more conservative choice, the latter (default) seems more
  accurate, but may underestimate the variance.}
}
\value{
  A data.frame with columns Expected and Variance
}
\description{
  Predict risk with missing data
}
\author{
  mg14
}

