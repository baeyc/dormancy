#' Chilling units
#'
#' Compute chilling units accumulated by a plant, according to a beta-shaped function.
#'
#' @name chillingUnits
#' @aliases chilling CU
#'
#' @param temp.data a data frame with the dates and temperatures to accumulate
#' @param var.names the name of the temperature and date variables, in the format \code{list(temp="temp.name",date="date.name",duration="duration.name")}
#' @param temp.min the threshold above which the plant is accumulating chilling units
#' @param temp.max the threshold below which the plant is accumulating chilling units
#' @param mu the temperature at which the plant accumulates the highest number of chilling units
#' @param s the sample size of the underlying Beta distribution
#' @return a vector with the chilling units corresponding to each input temperature
#'
#' @details The accumulation of chilling units is expressed as the Beta law probability distribution function (pdf). The underlying Beta law is re-parametrized
#' using the mode and the so-called sample size. More precisely, using the following definition for the pdf of the Beta law with parameters \code{a} and \code{b}:
#' \deqn{\frac{\Gamma(a+b)}{\Gamma (a)+\Gamma (b)} \ x^{a-1}(1-x)^{1-b}}
#' for \eqn{0 \leq x \leq 1}, \eqn{a>0} and \eqn{b >0}. Now, we further assume here that \eqn{a > 1} and \eqn{b > 1}, so that the mode exists, and we parametrized
#' the Beta distribution using \eqn{\mu} and \eqn{s} defined as:
#' \deqn{\mu = \frac{a-1}{a+b-2}, \quad s = a+b}
#'
#' Since the Beta law has its support on \eqn{[0,1]}, it is evaluated at \code{(x-temp.min)/(temp.max-temp.min)}, where \code{x} is the recorded temperature, and
#' \code{mu} is transformed into \code{(mu-temp.min)/(temp.max-temp.min)} to lie between 0 and 1
#'
#' @export chillingUnits
chillingUnits <- function(temp.data,
                          var.names=list(temp="temp",date="date",duration="duration"),
                          temp.min=-5,temp.max=15,mu,s){
  # go back to Beta law parameters
  mu <- (mu-temp.min)/(temp.max-temp.min) # mu is given in the [temp.min,temp.max] scale
  a <- mu*(s-2) + 1
  b <- (s-2)*(1-mu) + 1

  x <- (temp.data[,var.names$temp]-temp.min)/(temp.max-temp.min) # to lies between 0 and 1

  cu <- temp.data[,var.names$duration] * dbeta(x,shape1 = a, shape2 = b)

  return(cu)
}

