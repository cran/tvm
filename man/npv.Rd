\name{npv}
\alias{npv}
\title{Net Present Value of a cashflow (NPV)}
\usage{
npv(i, cf, ts = seq(from = 0, by = 1, along.with = cf))
}
\arguments{
  \item{i}{The rate used to discount the cashflow. It must
  be effective and with a periodicity that matches that of
  the cashflow}

  \item{cf}{The cashflow}

  \item{ts}{The times on which the cashflow ocurrs. It is
  assumed that \code{cf[idx]} happens at moment
  \code{t[idx]}}
}
\description{
Net Present Value of a cashflow (NPV)
}
\examples{
npv(i = 0.01, cf = c(-1, 0.5, 0.9), ts = c(0, 1, 3))
}

