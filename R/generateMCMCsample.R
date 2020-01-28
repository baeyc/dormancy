#' Generate a MCMC sample
generateMCMC<-function(current.state,obs.data,mcmc.params)
{
  p <- length(current.state)

  next.state <- rep(0,p)
  candidate <- rep(0,p)

  sd <- sqrt(mcmc.params$lambda*mcmc.params$var.proposal)

  if(mcmc.params$proposal == "GlRW")
  {
    candidate <- mvtnorm::rmvnorm(1, mcmc.params$mean.proposal, mcmc.params$var.proposal)

    # MH ratio
    joint.dist.candidate <- prod(modelProbaBB(temp.data = temp.data,
                                              var.names = var.names,
                                              origin.date = origin.date,
                                              temp.params = temp.params,
                                              stats::setNames(as.list(chains[[i]][1,]),names(theta.current)))$probaBB)


    r = min(1,ratioRW(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))
    probaAcc<-r
    #print(probaAcc)

    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    } else {
      next_state<-last_state
    }
  }
  else if(algo == "hgs")
  {
    if(prop == "CW")
    {
      candidate = unlist(last_state[1,])
      # we randomly choose one of the components
      p = sample(1:(size_state-1),1)
      if(p <= 2){
        if (dist == "unif"){
          a<-max(0,unlist(last_state[p]-sqrt(3)*sd[p]))
          b<-min(1,unlist(last_state[p]+sqrt(3)*sd[p]))
          candidate[p] <- runif(1,min=a,max=b)
        }else if(dist == "norm"){
          candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[p])
        }
      }else if(p < (size_state-1)){
        if (dist == "unif"){
          a<-last_state[p]-sqrt(3)*sd[p]
          b<-last_state[p]+sqrt(3)*sd[p]
          candidate[p] <- runif(1,min=a[[1]],max=b[[1]])
        }else if(dist == "norm"){
          candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[p])
        }
      }else{
        tauLast<-c(last_state$tauW,last_state$tauB)
        if (dist == "unif"){
          tau<-jointPropTauUnif(tauLast,sd[(size_state-1):size_state])
          candidate[(size_state-1):size_state]<-tau
        }else if(dist == "norm"){
          tau<-jointPropTauNorm(tauLast,sd[(size_state-1):size_state])
          candidate[(size_state-1):size_state]<-tau
        }
      }
      #print(paste("candidate :",candidate))
    }
    else if(prop == "Gl") ## check the differences with mh sampling ... mh should be multivariate?
    {
      a<-last_state-sqrt(3)*sd
      b<-last_state+sqrt(3)*sd
      a[[1]]<-max(0,a[[1]]);a[[2]]<-max(0,a[[2]]); # troncate mode to be between 0 and 1
      b[[1]]<-min(1,b[[1]]);b[[2]]<-min(1,b[[2]]); # troncate mode to be between 0 and 1
      for (p in 1:(size_state-1))
      {
        if (dist == "unif"){
          candidate[p] <- runif(1,min=a[[p]],max=b[[p]])
        }else if (dist == "norm")
          if(p <= 2){
            candidate[p] <- rtruncnorm(1,a=0,b=1,mean=last_state[[p]],sd=sd[[p]])
          }else{
            candidate[p] <- rnorm(1,mean=last_state[[p]],sd=sd[[p]])
          }
      }
      tauLast<-c(last_state$tauW,last_state$tauB)
      if (dist == "unif"){
        tau<-jointPropTauUnif(tauLast,sd[(size_state-1):size_state])
        candidate[(size_state-1):size_state]<-tau
      }else if(dist == "norm"){
        tau<-jointPropTauNorm(tauLast,sd[(size_state-1):size_state])
        candidate[(size_state-1):size_state]<-tau
      }
    }

    r = min(1,ratioRW(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))

    u<-runif(1)
    if(is.finite(r) & u <= r) {
      next_state<-candidate
    }else{
      next_state<-last_state
    }

    probaAcc<-rep(0,size_state)
    if(prop == "CW")
    {
      probaAcc[p]<-r
    }
    else if(prop == "Gl")
    {
      temp<-last_state
      for (j in 1:size_state)
      {
        temp[j] = candidate[j]
        r = min(1,ratioRW(candidate,last_state,data,hyper,copFamily,dist,sd[(size_state-1):size_state]))

        probaAcc[j] = r
        temp[j] = last_state[j]
      }
    }
  }
  #print(paste("r =",r,"probaAcc =",probaAcc))
  return(list(chain=next_state,lambda=lambda,accRate=probaAcc))
}
