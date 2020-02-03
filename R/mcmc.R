#' MCMC algorithm for Bayesian parameter estimation
#'
#' HELP PAGE TO WRITE --
#'
#' @name mcmc
#' @aliases mcmc
#'
#' @param data a list containing the plant and temperatures data
#' @param temp.params a list with the fixed parameters giving the minimum and maximum temperatures for computing chilling and forcing units
#' @param priors a list with the priors on the parameters
#' @param origin.date the date to be used as the origin when computing cumulative sum of temperatures
#' @param control a list of options for the algorithm
#'
#' @export mcmc
mcmc <- function(data=list(obs.data=obs.data, temp.data=temp.plants, var.names=var.names),
                   temp.params = list(temp.min.cu = -10, temp.max.cu = 15, temp.min.fu = 5, temp.max.fu = 35),
                   priors = list(a.cu = prior(distRNG="runif", hyperParams=list(min=-5, max=5)),
                                 b.cu = prior(distRNG="runif", hyperParams=list(min=2, max=200)),
                                 a.fu = prior(distRNG="runif", hyperParams=list(min=5, max=20)),
                                 b.fu = prior(distRNG="runif", hyperParams=list(min=1, max=20)),
                                 mu = prior(distRNG="rnorm", hyperParams=list(mean=1500, sd=1000)),
                                 s = prior(distRNG="rnorm", hyperParams=list(mean=750, sd=500))),
                   origin.date="09-01",
                   control = list(proposal="AdGl",
                                  size=100000)){
  # Initialize algorithm
  cat("Initialize MCMC algorithm ...")
  names.params <- names(priors)
  init.state <- sapply(1:length(names.params), FUN = function(i){do.call(priors[[i]]@distRNG, c(list(n = 1), priors[[i]]@hyperParams))})
  names(init.state) <- names.params
  p <- length(names.params)
  print(init.state)

  mean.chain <- init.state
  var.chain <- diag((2.38/sqrt(p))*rep(1,p))
  alpha.optim <- ifelse(control$proposal=="AdGl", 0.234, 0.44)
  lambda <- rep(1,p)
  stoch.step.power <- 0.6
  accept.rates <- matrix(rep(0,p),nr=1,nc=p)

  # Compute likelihood and priors values for initial state
  pij <- modelProbaBB(temp.data = data$temp.data,
                      var.names = data$var.names,
                      origin.date = origin.date,
                      temp.params = temp.params,
                      stats::setNames(as.list(init.state),names.params))
  pij <- dplyr::inner_join(data$obs.data,pij,by = c("session", "plant", "rep"))
  likelihood.current <- exp(sum(dbinom(pij$budburst,1,pij$probaBB,log = TRUE)))
  prior.current <- sapply(1:length(names.params), FUN = function(i){dname <- priors[[i]]@distRNG;
                                                                    substr(dname,1,1) <- "d"
                                                                    do.call(dname, c(list(x = init.state[i]), priors[[i]]@hyperParams))})
  prior.current <- prod(prior.current)

  # Generate MCMC
  chain <- matrix(init.state,nrow=1,ncol=p)
  pb <- txtProgressBar(min=1,max=control$size,style = 3)
  for (m in 1:control$size){
    setTxtProgressBar(pb,m)
    # generate candidate
    if (control$proposal == "AdGl"){
      candidate <- as.vector(mvtnorm::rmvnorm(1, chain[m,], diag(lambda)%*%var.chain)) # in the AdGl case, lambda has the same value on all its components
    }else if (control$proposal == "AdCW"){
      lambda.sqrt.mat <- diag(sqrt(lambda))
      candidate <- as.vector(mvtnorm::rmvnorm(1, chain[m,], lambda.sqrt.mat%*%var.chains%*%lambda.sqrt.mat))
    }

    # MH ratio
    pij <- modelProbaBB(temp.data = data$temp.data,
                        var.names = data$var.names,
                        origin.date = origin.date,
                        temp.params = temp.params,
                        stats::setNames(as.list(candidate),names.params))
    pij <- dplyr::inner_join(data$obs.data,pij,by = c("session", "plant", "rep"))
    likelihood.candidate <- exp(sum(dbinom(pij$budburst,1,pij$probaBB,log = TRUE)))
    prior.candidate <- sapply(1:length(names.params), FUN = function(i){dname <- priors[[i]]@distRNG;
                                                                        substr(dname,1,1) <- "d"
                                                                        do.call(dname, c(list(x = candidate[i]), priors[[i]]@hyperParams))})
    prior.candidate <- prod(prior.candidate)

    ratio <- (likelihood.candidate*prior.candidate)/(likelihood.current*prior.current) # symetric RW

    # Set next state of the chain
    if (is.na(ratio)){
      next.state <- chain[m,]
    }else{
      u <- runif(1)
      if (u<=ratio){
        next.state <- candidate
        likelihood.current <- likelihood.candidate
        prior.current <- prior.candidate
      }else{
        next.state <- chain[m,]
      }
    }

    accept.rates <- rbind(accept.rates,(accept.rates[m,]*(m-1) + (!is.na(ratio) & u <= ratio))/(m))
    chain <- rbind(chain,next.state)

    # Update params of MCMC sampler
    stoch.step <- 1/(m^stoch.step.power)
    mean.chain <- mean.chain + stoch.step * (next.state - mean.chain)
    var.chain <- var.chain + stoch.step * ((next.state-mean.chain)%*%t(next.state-mean.chain) - var.chain)

    lambda <- lambda * exp(stoch.step*(accept.rates[m+1,] - alpha.optim))
  }

  return(list(chain=chain,ar=accept.rates,lambda=lambda))
}
