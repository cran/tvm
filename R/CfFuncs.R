#' @useDynLib tvm

#' @title Adjusts the discount factors by a spread
#' 
#' @param fd vector of discount factors used to discount cashflows in \code{1:length(fd)} periods
#' @param spread effective spread
#' @export
adjust_disc <- function(fd,spread) {
  zeros <- (1/fd)^(1/seq(along.with=fd))
  zeros_adj <- zeros + spread
  1/(zeros_adj^(seq(along.with=zeros_adj)))
}

#' Calculates the Total Financial Cost (CFT)
#' 
#' This is the IRR of the loan's cashflow, after adding all the extra costs
#' 
#' It is assumed that the loan has monthly payments
#' The CFT is returned as an effective rate of periodicty equal to that of the maturity and the rate
#' The interest is calculated over amt + fee
#' 
#' @param amt The amount of the loan
#' @param maturity The maturity of the loan
#' @param rate The loan rate, in  effective rate
#' @param up_fee The fee that the loan taker pays upfront
#' @param per_fee The fee that the loan payer pays every period
#' @export
cft <- function(amt, maturity, rate, up_fee = 0, per_fee = 0) {
  full_amt <- amt + up_fee
  p <- pmt(amt = full_amt, maturity = maturity, rate = rate)
  rate(amt = amt, maturity = maturity, pmt = p + per_fee)
}

#' Net Present Value of a cashflow (NPV)
#' 
#' @param i The rate used to discount the cashflow. It must be effective and with a periodicity that matches that of the cashflow
#' @param cf The cashflow
#' @param t The times on which the cashflow ocurrs. It is assumed that \code{cf[idx]} happens at moment \code{t[idx]}
#' @export
npv <- function(i, cf, t=seq(from=0,by=1,along.with=cf)) sum(cf/(1+i)^t)

#' Internal Rate of Return of a cashflow (IRR)
#' 
#' The IRR is returned as an effective rate with periodicity equal to that of the cashflow
#' 
#' @param cf The cashflow
#' @param t The times on which the cashflow ocurrs. It is assumed that \code{cf[idx]} happens at moment \code{t[idx]}
#' @export
irr <- function(cf, t=seq(from=0,by=1,along.with=cf)) { uniroot(npv, c(0,100000), cf=cf, t=t)$root }

#' The value of the payment of a loan with constant payments (french type amortization)
#' 
#' The periodicity of the maturity and the rate must match, and this will be the periodicity of the payments
#'
#' @param amt The amount of the loan
#' @param maturity The maturity of the loan
#' @param rate The rate of the loan
#' @export
pmt <- function(amt, maturity, rate) {  
  return(amt*rate/(1-(1+rate)^(-maturity)))
}

#' The rate of a loan with constant payments (french type amortization)
#' 
#' The periodicity of the maturity and the payment must match, and this will be the periodicity of the rate (which is returned as an effective rate)
#' 
#' @param amt The amount of the loan
#' @param maturity The maturity of the loan
#' @param pmt The payments of the loan
#' @param extrema Vector of length 2 that has the minimum and maximum value to search for the rate
#' @param tol The tolerance to use in the root finding algorithm
#' @export
rate <- function(amt, maturity, pmt, extrema=c(1e-4,1e9), tol=1e-4) {   
  zerome <- function(r) amt/pmt-(1-1/(1+r)^maturity)/r
  if(zerome(extrema[1])>0) return(0)
  if(zerome(extrema[2])<0) return(extrema[2])
  return(uniroot(zerome, interval=extrema, tol=tol)$root)
}

rate <- Vectorize(FUN=rate,vectorize.args=c("amt","maturity","pmt"))

#' Constructor for the loan class
#' 
#' Creates a loan
#' 
#' @param rate The periodic effective rate of the loan
#' @param maturity The maturity of the loan, measured in the same units as the periodicity of the rate
#' @param amt The amount loaned
#' @param type The type of loan. Available types are \code{c("bullet","french","german")}
#' @param grace_int The number of periods that the loan doesn't pay interest and capitalizes it. Leave in 0 for zero loans
#' @param grace_amort The number of periods that the loan doesn't amortize
#' @export
loan <- function(rate,maturity,amt,type,grace_int = 0, grace_amort = grace_int) {
  stopifnot(grace_int <= grace_amort)
  stopifnot(grace_amort < maturity)
  l <- structure(list(rate = rate, maturity = maturity, amt = amt, type = type, grace_amort = 0, grace_int = 0), class = c(type,"loan"))
  if (grace_amort > 0 || grace_int > 0) {
    sl <- loan(rate=rate,maturity=maturity-grace_amort,amt=1,type=type)
    l$cf <- amt*(1+rate)^grace_int*c(rep_len(0,grace_int),rep_len(rate,grace_amort-grace_int),sl$cf)  
  } else {
    l$cf <- cashflow(l)  
  }
  l
}

cashflow <- function(l) {
  UseMethod("cashflow")
}

cashflow.loan <- function(l) {
  stop("Can't get cashflow for a loan without the proper type")
}

#' Cashflow for a bullet loan
#' 
#' A bullet loan pays interest periodically and returns its principal at maturity
#' 
#' @param l The loan
cashflow.bullet <- function(l) {
  f <- rep_len(l$rate,l$maturity)
  f[l$maturity] <- 1 + f[l$maturity]
  f*l$amt
}

#' Cashflow for a german loan
#' 
#' A german loan has constant amortizations
#' 
#' @param l The loan
cashflow.german <- function(l) {  
  k <- rep_len(1/l$maturity,l$maturity)
  krem <- c(1,1-cumsum(k))
  i <- head(krem,l$maturity)*l$rate
  (k+i)*l$amt
}

#' Cashflow for a french loan
#' 
#' A french loan has a constant payment
#' 
#' @param l The loan
cashflow.french <- function(l) {  
  rep_len(l$rate / (1 - (1+l$rate)^(-l$maturity)),l$maturity)*l$amt 
}

#' Value of a discounted cashflow
#' 
#' @param disc The discount factor curve
#' @param cf The cashflow
#' @export
disc_cf <- function(disc,cf) {
  sum(disc*cf)
}

#' Remaining capital in a loan
#' 
#' The amount that has to be repayed at each moment in a loan, at the end of the period
#' 
#' @param cf The cashflow of the loan
#' @param amt The original amount of the loan
#' @param r The periodic rate of the loan
#' @export
rem <- function(cf,amt,r) {
  s <- function(t) amt*(1+r)^t-sum(cf[1:t]*(1+r)^(t-(1:t)))
  vapply(X=seq_along(cf),FUN=s,FUN.VALUE=1)
}