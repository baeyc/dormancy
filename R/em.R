#' EM algorithm for probability of budburst
#'
#' HELP PAGE TO WRITE -- EM algorithme for estimation of the parameters when they vary accross years
#'
#' @name algoEM
#' @aliases algoEM
#'
#' @param data a list containing the plant and temperatures data
#' @param temp.params a list with the fixed parameters giving the minimum and maximum temperatures for computing chilling and forcing units
#' @param init.params a list with the initial values for the parameters to be estimated
#' @param origin.date the date to be used as the origin when computing cumulative sum of temperatures
#' @param control a list of options for the algorithm
#'
#' @export algoEM
algoEM <- function(data,
                   temp.params = list(temp.min.cu = -10, temp.max.cu = 15, temp.min.fu = 5, temp.max.fu = 35),
                   init.params = list(a.cu = -2, b.cu = 5, a.fu = 15, b.fu = 4, mu = 1500, s = 750),
                   origin.date="09-01",
                   control = list(proposal="AdGl",
                                  mcsize=5,
                                  saem.size1=10,
                                  saem.size2=0)){
  # Initialize algorithm
  beta.current <- unlist(init.params)
  names.params <- names(beta.current)
  Gamma.current <- diag(length(beta.current))

  # Initialize Markov Chains and MCMC sampler (one per session)
  n <- length(unique(obs.data$session))
  p <- length(init.params)
  chains <- lapply(1:n, FUN = function(i){tmp<-matrix(beta.current,nrow=1,ncol=p);colnames(tmp)<-names(beta.current);return(tmp)})
  mean.chains <- lapply(1:n,FUN = function(i){beta.current})
  var.chains <- lapply(1:n,FUN = function(i){diag((2.38/sqrt(p))*rep(1,p))})
  alpha.optim <- ifelse(control$proposal=="AdGl", 0.234, 0.44)
  mcmc.params <- list(mean.chains=mean.chains,
                      var.chains=var.chains,
                      lambda=lapply(1:n,FUN=function(i){rep(1,p)}),
                      stoch.step.power=0.6,
                      accept.rates=lapply(1:n,FUN=function(i){matrix(rep(0,p),nr=1,nc=p)}),
                      proposal=control$proposal)

  cond.dist.current <- sapply(1:n,FUN = function(i){pij <- modelProbaBB(temp.data = data$temp.data,
                                                                       var.names = data$var.names,
                                                                       origin.date = origin.date,
                                                                       temp.params = temp.params,
                                                                       stats::setNames(as.list(chains[[i]][1,]),names.params));
                                                     pij <- dplyr::inner_join(obs.data,pij);
                                                     exp(sum(pij$budburst*log(pij$probaBB)) + sum((1-pij$budburst)*log(1-pij$probaBB)))
                                                     })

  # Initialize sufficient statistics
  t1 <- lapply(1:n, FUN = function(i){ (-0.5)*t(chains[[i]])%*%chains[[i]]}) # for Gamma^-1 beta
  t1 <- Reduce('+',t1)
  t2 <- Reduce('+',chains)

  # Step in the SAEM algorithm (to be added to the control option)
  stoch.step.SAEM <- 1
  alpha <- 0.6

  # SAEM loop
  for (k in 1:(control$saem.size1+control$saem.size2)){

    # (SA) E-step

    # Generate MCMC
    #no.cores <- max(1,parallel::detectCores() - 1)
    # Initiate cluster
    #doParallel::registerDoParallel(no.cores)

    #per.indiv <- foreach::foreach(i=1:n) %dopar% {
    per.indiv <- list()
    for (i in 1:n){
      chain.indiv <- tail(chains[[i]],1)
      ar.indiv <- tail(mcmc.params$accept.rates[[i]],1)
      mean.chains <- mcmc.params$mean.chains[[i]]
      var.chains <- mcmc.params$var.chains[[i]]
      lambda <- mcmc.params$lambda[[i]]
      for (m in 1:control$mcsize){
        # generate candidate
        if (control$proposal == "AdGl"){
          candidate <- mvtnorm::rmvnorm(1, tail(chain.indiv,1), diag(lambda)%*%var.chains) # in the AdGl case, lambda has the same value on all its components
        }else if (control$proposal == "AdCW"){
          lambda.sqrt.mat <- diag(sqrt(lambda))
          candidate <- mvtnorm::rmvnorm(1, tail(chain.indiv,1), lambda.sqrt.mat%*%var.chains%*%lambda.sqrt.mat)
        }

        # MH ratio
        pij <- modelProbaBB(temp.data = temp.data,
                            var.names = var.names,
                            origin.date = origin.date,
                            temp.params = temp.params,
                            stats::setNames(as.list(candidate),names.params));
        pij <- dplyr::inner_join(obs.data,pij);
        cond.dist.candidate <- exp(sum(pij$budburst*log(pij$probaBB)) + sum((1-pij$budburst)*log(1-pij$probaBB)))

        ratio <- (cond.dist.candidate/cond.dist.current[i]) * mvtnorm::dmvnorm(candidate,tail(chains[[i]],1),Gamma.current)

        # Set next state of the chain
        if (is.na(ratio)) next.state <- tail(chain.indiv,1)
        u <- runif(1)
        if (!is.na(ratio)) next.state <- (u<=ratio)*candidate + (u>ratio)*tail(chain.indiv,1)
        iter <- (k-1)*control$mcsize+m
        ar.indiv <- rbind(ar.indiv,(tail(ar.indiv,1)*(iter-1) + (!is.na(ratio) & u <= ratio))/(iter))
        chain.indiv <- rbind(chain.indiv,next.state)

        # Update params of MCMC sampler
        mean.chains <- mean.chains + mcmc.params$stoch.step * (next.state - mean.chains)
        var.chains <- var.chains + mcmc.params$stoch.step *
          (t(next.state-mean.chains)%*%(next.state-mean.chains) - var.chains)

        lambda <- as.vector(lambda * exp((1/((k-1)*control$mcsize+m)^mcmc.params$stoch.step.power)*(tail(ar.indiv,1) - alpha.optim)))
      }

      #list(mean.chains=mean.chains,var.chains=var.chains,lambda=lambda,accept.rates=ar.indiv,chain.indiv=chain.indiv)
      per.indiv[[i]] <- list(mean.chains=mean.chains,var.chains=var.chains,lambda=lambda,accept.rates=ar.indiv,chain.indiv=chain.indiv)
    }

    # Store the results in the original objects
    chains <- lapply(1:n, FUN=function(i){rbind(chains[[i]],tail(per.indiv[[i]]$chain.indiv,control$mcsize))})
    mcmc.params$accept.rates <- lapply(1:n,FUN=function(i){rbind(mcmc.params$accept.rates[[i]],tail(per.indiv[[i]]$accept.rates,5))})
    mcmc.params$mean.chains <- lapply(1:n,FUN=function(i){per.indiv[[i]]$mean.chains})
    mcmc.params$var.chains <- lapply(1:n,FUN=function(i){per.indiv[[i]]$var.chains})
    mcmc.params$lambda <- lapply(1:n,FUN=function(i){per.indiv[[i]]$lambda})

    # update sufficient statistics
    t1.new <- lapply(1:n, FUN = function(i){
      averageMCMC <- lapply(1:control$mcsize, FUN=function(j){
        (-0.5)*chains[[i]][control$mcsize-j+1,]%*%t(chains[[i]][control$mcsize-j+1,])})
      return(Reduce('+',averageMCMC)/control$mcsize)
      })
    t1.new <- Reduce('+',t1.new)
    t2.new <- lapply(1:n, FUN = function(i){
      averageMCMC <- lapply(1:control$mcsize, FUN=function(j){chains[[i]][control$mcsize-j+1,]})
      return(Reduce('+',averageMCMC)/control$mcsize)
    })
    t2.new <- Reduce('+',t2.new)

    # Stochastic approximation
    t1 <- t1 + stoch.step.SAEM*(t1.new-t1)
    t2 <- t2 + stoch.step.SAEM*(t2.new-t2)

    # M-step
    beta.current <- t2/n
    Gamma.current <- t1/n - t(beta.current)%*%beta.current

    # Stochastic step
    if (k>control$saem.size1){stoch.step.SAEM <- (1/(k+control$saem.size1))^alpha}
  }
}
