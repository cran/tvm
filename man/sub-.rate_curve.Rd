\name{[.rate_curve}
\alias{[.rate_curve}
\title{Returns a particular rate or rates from a curve}
\usage{
\method{[}{rate_curve}(r, rate_type = "zero", x = NULL)
}
\arguments{
  \item{r}{The rate_curve object}

  \item{x}{The points in time to return}

  \item{rate_type}{The rate type}
}
\value{
If \code{x} is \code{NULL}, then returns a rate function of
\code{rate_type} type. Else, it returns the rates of
\code{rate_type} type and corresponding to time \code{x}
}
\description{
Returns a particular rate or rates from a curve
}
\examples{
r <- rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero")
r["zero"]
r["swap",c(1.5, 2)]
}

