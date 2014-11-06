\name{MakeTimeDependent}
\alias{MakeTimeDependent}
\title{Split a data.frame by a time-dependent variable for use in a coxph model}
\usage{
  MakeTimeDependent(dataFrame, timeEvent,
    timeStop = dataFrame$time,
    timeStart = rep(0, nrow(dataFrame)),
    status = dataFrame$status)
}
\arguments{
  \item{dataFrame}{A data frame}

  \item{timeEvent}{The time of the time-deependent
  covariate (assumed to be a simple event, such as a
  transplant).}

  \item{timeStop}{Start time for the survival object.
  Default: rep(0, nrow(dataFrame)).}

  \item{timeStart}{Stop time of the survival object.
  Default: dataFrame$time.}

  \item{status}{Status of the survival object (0=alive,
  1=dead). Default: dataFrame$status.}
}
\value{
  A data.frame with extra extra rows for observations after
  the event and extra columns time1, time2 (start and stop
  times), event for the event.
}
\description{
  Split a data.frame by a time-dependent variable for use
  in a coxph model
}
\author{
  mg14
}

