\name{SimSurvNonp}
\alias{SimSurvNonp}
\title{Non-parametric survival simulations}
\usage{
  SimSurvNonp(risk, surv, H0 = basehaz(coxph(surv ~ 1)))
}
\arguments{
  \item{risk}{A vector of the relative risk/}

  \item{surv}{A survival object}

  \item{H0}{The baseline hazed function. By default this is
  given by basehaz(coxph(surv ~ 1)).}
}
\value{
  A Surv()ival object
}
\description{
  Non-parametric survival simulations
}
\author{
  mg14
}

