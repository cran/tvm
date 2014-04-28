\name{irr}
\alias{irr}
\title{The IRR is returned as an effective rate with periodicity equal to that of the cashflow}
\usage{
irr(cf, ts = seq(from = 0, by = 1, along.with = cf), interval = c(0, 1e+05))
}
\arguments{
  \item{cf}{The cashflow}

  \item{ts}{The times on which the cashflow ocurrs. It is
  assumed that \code{cf[idx]} happens at moment
  \code{t[idx]}}

  \item{interval}{A length 2 vector that indicates the root
  finding algorithm where to search for the irr}
}
\description{
Internal Rate of Return of a cashflow (IRR)
}
\examples{
irr(cf = c(-1, 0.5, 0.9), ts = c(0, 1, 3))
}

