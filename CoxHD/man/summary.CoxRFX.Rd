\name{summary.CoxRFX}
\alias{summary.CoxRFX}
\title{A summary method for CoxRFX models}
\usage{
  summary.CoxRFX(x)
}
\arguments{
  \item{x}{A CoxRFX model}
}
\value{
  NULL
}
\description{
  This model prints the means and variances for each groups
  of covariates, as well as the variance components. For
  the means a Wald test with 1df is computed testing the
  null-hypothesis of being zero. The null-hypothesis of
  zero variance is tested using a combined Wald test that
  all coefficients in the group are identical to the mean.
  This test has tr(H-1 I) df.
}
\author{
  mg14
}

