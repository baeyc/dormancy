#' Class "prior"
#'
#' An object of the \code{prior} class
#'
#' @name prior-class
#' @docType class
#' @aliases prior-class prior
#' @slot distName the name of the prior distribution
#' @slot hyperParams a list of the distribution parameters (hyperparameters)

prior <- setClass("prior",
                  slots = c(distRNG="character",
                            hyperParams="list"
                  ),

                  # Set the default values for the slots
                  prototype=list(
                    distRNG = "fixed",
                    hyperParams = list()
                  ),

                  validity=function(object)
                  {
                    if (!object@distRNG%in%c("fixed","rnorm","rbeta","rlnorm","rlnormTrunc","rchisq","rexp","rf","rgamma","rt","runif","rweibull","rtruncnorm")){
                      stop("Chosen prior is not currently supported. Please choose between: rnorm, rbeta, rlnorm, rchisq, rexp, rf, rgamma, rt, runif, rweibull, rtruncnorm")
                    }
                    return(TRUE)
                  }
)

