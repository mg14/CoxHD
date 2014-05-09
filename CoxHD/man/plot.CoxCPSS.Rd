\name{plot.CoxCPSS}
\alias{plot.CoxCPSS}
\title{Plot a CoxCPSS model}
\usage{
  plot.CoxCPSS(stabCox, xlab = "1/lambda",
    ylab = "Selection probability",
    lty = rep(1, ncol(stabCox$X)), log = "",
    col = c("black", "grey"),
    xlim = range(1/stabCox$lambda[stabCox$Lambda]),
    ylim = c(0, 1), ...)
}
\arguments{
  \item{stabCox}{}

  \item{xlab}{}

  \item{ylab}{}

  \item{lty}{}

  \item{log}{}

  \item{col}{}

  \item{xlim}{}

  \item{ylim}{}

  \item{...}{}
}
\value{
  NULL
}
\description{
  Plot a CoxCPSS model
}
\author{
  mg14
}

