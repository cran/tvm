\name{plot.rate_curve}
\alias{plot.rate_curve}
\title{Plots a rate curve}
\usage{
\method{plot}{rate_curve}(x, rate_type = NULL, ...)
}
\arguments{
  \item{x}{The rate curve}

  \item{rate_type}{The rate types to plot, in
  c("french","fut","german","zero","swap")}

  \item{...}{Other arguments (unused)}
}
\description{
Plots a rate curve
}
\examples{
r <- rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero")
plot(r)
\dontrun{
plot(r, rate_type = "german")
plot(r, rate_type = c("french", "german"))
}
}

