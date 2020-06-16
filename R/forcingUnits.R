#' Forcing units
#'
#' Compute forcing units accumulated by a plant, according to a beta-shaped function.
#'
#' @name forcingUnits
#' @aliases forcing CU
#'
#' @param temp.data a data frame with the dates and temperatures to accumulate
#' @param var.names the name of the temperature and date variables, in the format \code{list(temp="temp.name",date="date.name",duration="duration.name")}
#' @param temp.min the threshold above which the plant is accumulating chilling units
#' @param temp.max the threshold below which the plant is accumulating chilling units
#' @param a the temperature at which the plant accumulates half the maximum of forcing units
#' @param b the rate of forcing units accumulation
#' @return a vector with the forcing units corresponding to each input temperature
#'
#' @details The accumulation of forcing units is expressed through a logistic curve with parameters \code{a} and \code{b}:
#' \deqn{fu(t) = \frac{1}{1 + \exp(-(t-a)/b)} }
#' if \eqn{temp.min \leq t \leq temp.max}, and 0 otherwise, and where \eqn{a} is the temperature at which the plant accumulates half the maximum of forcing units and \eqn{b} is the rate of accumulation
#'
#' @export forcingUnits
forcingUnits <- function(temp.data,
                         var.names=list(temp="temp",date="date",duration="duration"),
                         temp.min=5,temp.max=30,a,b){

  x <- unlist(temp.data[,var.names$temp])
  fu <- temp.data[,var.names$duration] * ifelse(temp.min <= x & x <= temp.max,1/(1+exp(-(x-a)/b)),0)

  return(unlist(fu))

}
