\name{WaldTest}
\alias{WaldTest}
\title{A Wald test for the coefficients of a CoxRFX model}
\usage{
  WaldTest(coxRFX, var = c("var2", "var"))
}
\arguments{
  \item{coxRFX}{The CoxRFX model}

  \item{var}{Which type of variance estimate to use. The
  default choice is var2 = H-1 I H-1. A more conservative
  choice is var2 = H-1.}
}
\value{
  A data.frame with columns coef, sd, z and p
}
\description{
  This separately tests the null-hypothesis of being zero
  on each coefficient in a CoxRFX model using a Wald test.
  The test statistic is z^2 = coef^2/var[coef].
}
\author{
  mg14
}

