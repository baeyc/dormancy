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
                   control = list(proposal="AdGl",
                                  size=100000),
                   continue = FALSE,
                   last.state = NULL){
  # Initialize algorithm
  cat("Initialize MCMC algorithm...\n")
  names.params <- names(priors)
  init.state <- sapply(1:length(names.params), FUN = function(i){do.call(priors[[i]]@distRNG, c(list(n = 1), priors[[i]]@hyperParams))})
  names(init.state) <- names.params
  p <- length(names.params)
  print(init.state)

  proposal <- control$proposal

  mean.chain <- init.state
  var.chain <- diag((2.38/sqrt(p))*rep(1,p))
  lambda <- rep(1,p)
  stoch.step.power <- 0.6
  accept.rates <- matrix(rep(0,p),nr=1,nc=p)

  # Compute likelihood and priors values for initial state
  pij <- modelProbaBB(temp.data = data$temp.data,
                      var.names = data$var.names,
                      temp.params = temp.params,
                      cufu.params = stats::setNames(as.list(init.state),names.params))
  pij <- dplyr::inner_join(data$obs.data,pij,by = c("session", "plant", "rep"))
  likelihood.current <- exp(sum(dbinom(pij$budburst,1,pij$probaBB,log = TRUE)))
  prior.current <- sapply(1:length(names.params), FUN = function(i){dname <- priors[[i]]@distRNG;
                                                                    substr(dname,1,1) <- "d"
                                                                    do.call(dname, c(list(x = init.state[i]), priors[[i]]@hyperParams))})
  prior.current <- prod(prior.current)


  # if we continue the chain from a previous stopped state
  if (continue){
    init.state <- tail(last.state$chain,1)
    mean.chain <- last.state$mean.and.var[[1]]
    var.chain <- last.state$mean.and.var[[2]]

    lambda <- last.state$lambda
    accept.rates <- tail(last.state$ar,1)

    # Compute likelihood and priors values for initial state
    pij <- modelProbaBB(temp.data = data$temp.data,
                        var.names = data$var.names,
                        temp.params = temp.params,
                        cufu.params =stats::setNames(as.list(init.state),names.params))
    pij <- dplyr::inner_join(data$obs.data,pij,by = c("session", "plant", "rep"))
    likelihood.current <- exp(sum(dbinom(pij$budburst,1,pij$probaBB,log = TRUE)))
    prior.current <- sapply(1:length(names.params), FUN = function(i){dname <- priors[[i]]@distRNG;
    substr(dname,1,1) <- "d"
    do.call(dname, c(list(x = init.state[i]), priors[[i]]@hyperParams))})
    prior.current <- prod(prior.current)
  }

  # Generate MCMC
  chain <- matrix(init.state,nrow=1,ncol=p)
  pb <- txtProgressBar(min=1,max=control$size,style = 3)
  for (m in 1:control$size){
    setTxtProgressBar(pb,m)

    if (control$proposal == "random") proposal <- sample(c("AdGl","CWAdCW","GlAdCW"),1)

    alpha.optim <- ifelse(proposal=="AdGl", 0.234, 0.44)

    # generate candidate
    ratio <- NA
    while(is.na(ratio) | is.infinite(ratio)){
      if (proposal == "AdGl"){
        candidate <- as.vector(mvtnorm::rmvnorm(1, chain[m,], lambda[1]*var.chain)) # in the AdGl case, lambda has the same value on all its components
      }else if (proposal == "GlAdCW"){
        lambda.sqrt.mat <- diag(sqrt(lambda))
        candidate <- as.vector(mvtnorm::rmvnorm(1, chain[m,], lambda.sqrt.mat%*%var.chain%*%lambda.sqrt.mat))
      }else if (proposal == "CWAdCW"){
        k <- sample(p,1)
        candidate <- chain[m,]
        candidate[k] <- candidate[k] + rnorm(1,0,sqrt(lambda[k]*var.chain[k,k]))
      }

      # MH ratio
      pij <- modelProbaBB(temp.data = data$temp.data,
                          var.names = data$var.names,
                          temp.params = temp.params,
                          stats::setNames(as.list(candidate),names.params))
      pij <- dplyr::inner_join(data$obs.data,pij,by = c("session", "plant", "rep"))
      likelihood.candidate <- exp(sum(dbinom(pij$budburst,1,pij$probaBB,log = TRUE)))
      prior.candidate <- sapply(1:length(names.params), FUN = function(i){dname <- priors[[i]]@distRNG;
                                                                          substr(dname,1,1) <- "d"
                                                                          do.call(dname, c(list(x = candidate[i]), priors[[i]]@hyperParams))})
      prior.candidate <- prod(prior.candidate)

      ratio <- min(1,(likelihood.candidate*prior.candidate)/(likelihood.current*prior.current)) # symetric RW
    }

    # Set next state of the chain
    u <- runif(1)
    if (u<=ratio){
      next.state <- candidate
      likelihood.current <- likelihood.candidate
      prior.current <- prior.candidate
    }else{
      next.state <- chain[m,]
    }

    accept.rates <- rbind(accept.rates,1*(next.state!=chain[m,]))
    chain <- rbind(chain,next.state)

    # Update params of MCMC sampler
    stoch.step <- 1/((m+1)^stoch.step.power)
    var.chain <- var.chain + stoch.step * ((next.state-mean.chain)%*%t(next.state-mean.chain) - var.chain)
    mean.chain <- mean.chain + stoch.step * (next.state - mean.chain)

    # Adapting lambda
    if (proposal == "AdGl"){
      lambda <- lambda * exp(stoch.step*(ratio - alpha.optim))
    }else if (proposal == "CWAdCW"){
      lambda[k] <- lambda[k] * exp(stoch.step*(ratio - alpha.optim))
    }else if (proposal == "GlAdCW"){
      for (k in 1:p){
        candidate.cw <- chain[m,]
        candidate.cw[k] <- candidate[k]
        pij <- modelProbaBB(temp.data = data$temp.data,
                            var.names = data$var.names,
                            temp.params = temp.params,
                            stats::setNames(as.list(candidate.cw),names.params))
        pij <- dplyr::inner_join(data$obs.data,pij,by = c("session", "plant", "rep"))
        likelihood.candidate.cw <- exp(sum(dbinom(pij$budburst,1,pij$probaBB,log = TRUE)))
        prior.candidate.cw <- sapply(1:length(names.params), FUN = function(i){dname <- priors[[i]]@distRNG;
                              substr(dname,1,1) <- "d"
                              do.call(dname, c(list(x = candidate.cw[i]), priors[[i]]@hyperParams))})
        prior.candidate.cw <- prod(prior.candidate.cw)

        ratio.cw <- min(1,(likelihood.candidate.cw*prior.candidate.cw)/(likelihood.current*prior.current)) # symetric RW

        lambda[k] <- lambda[k] * exp(stoch.step*(ratio.cw - alpha.optim))
      }
    }
  }

  return(list(chain=chain,ar=apply(accept.rates,2,cumsum)/nrow(accept.rates),lambda=lambda,mean.and.var=list(mean.chain,var.chain)))
}
