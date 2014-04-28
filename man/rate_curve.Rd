\name{rate_curve}
\alias{rate_curve}
\title{Creates a rate curve instance}
\usage{
rate_curve(rates = NULL, rate_type = "zero", pers = 1:length(rates),
  fun_d = NULL, fun_r = NULL, knots = 1:length(rates))
}
\arguments{
  \item{rates}{A rate vector}

  \item{rate_type}{The rate type. Must be on of
  c("fut","zero","swap")}

  \item{pers}{The periods the rates correspond to}

  \item{fun_d}{A discount factor function. fun_d(x) returns
  the discount factor for time x, vectorized on x}

  \item{fun_r}{A rate function. fun_r(x) returns the EPR
  for time x, vectorized on x}

  \item{knots}{The nodes used to bootstrap the rates}
}
\description{
Creates a rate curve instance
}
\note{
Currently a rate curve can only be built from one of the
following sources \enumerate{ \item A discount factor
function \item A rate function and a rate type from the
following types: "fut", "zero" or "swap" \item A rate
vector, a pers vector and a rate type as before }
}
\examples{
rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero")
rate_curve(fun_r = function(x) rep_len(0.1, length(x)), rate_type = "swap")
rate_curve(fun_d = function(x) 1 / (1 + x))
}

