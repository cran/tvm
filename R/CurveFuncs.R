#' Converts a futures curve to a spot curve
#' 
#' @param fut The futures curve
#' @export
fut_to_spot <- function(fut) {
  (cumprod(1+fut))^(1/(seq_along(fut)))-1
}

#' Converts a discount factor curve to a swap rate curve
#' 
#' @param disc The discount curve
#' @export
disc_to_swap <- function(disc) {
  (1-disc) / cumsum(disc)
}

#' Converts a swap curve to a discount factor curve
#' 
#' @param swap The swap curve
#' @export
swap_to_disc <- function(swap) {
  d = swap
  d[1] = 1 / (1 + swap[1])
  for (j in 2:length(swap)) {
    d[j] = (1 - sum(d[1:(j-1)])*swap[j]) / (1 + swap[j])
  }
  d
}

#' Converts a discount factor curve to a spot rate curve
#' 
#' @param disc The discount factor curve
#' @export
disc_to_spot <- function(disc) {
  (1 / disc)^(1/(seq_along(disc)))-1
}

#' Converts a spot curve to a futures curve
#' 
#' @param spot The spot curve
#' @export
spot_to_disc <- function(spot) {
  1 / ( (1 + spot)^(seq_along(spot)) )
}

#' Converts a discount factor curve to a futures curve
#' 
#' @param disc The discount factor curve
#' @export
disc_to_fut <- function(disc) {  
  exp(-diff(log(c(1,disc))))-1  
}

#' Converts a futures curve to a discount factor curve
#' 
#' @param fut The futures curve
#' @export
fut_to_disc <- function(fut) {
  1 / cumprod(1+fut)
}

#' Converts a discount factor curve to a german loan rate curve
#' 
#' @param disc The discount factor curve
#' @export
disc_to_german <- function(disc) {
  vapply(
    1:length(disc),
    function(i) ( 1 - 1 / i * sum(disc[1:i]) ) / sum( disc[1:i] * (1 - ( 1:i - 1 ) / i ) ),
    1)  
}

#' Converts a discount factor curve to a french loan rate curve
#' 
#' @param disc The discount factor curve
#' @export
disc_to_french <- function(disc) {  
  zerome = function(r,i,disc) 1/r * ( 1 - (1+r)^(-i) ) - sum(disc[1:i])
  vapply(
    1:length(disc),
    function(i) uniroot(function(r) zerome(r,i,disc),c(0.0001,1))$root,
    1)
}