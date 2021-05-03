# Analyze MCMC results

# Import data
dir <- "data/results/AdaptiveGlobal/Priors unif/"
files <- list.files(path = dir, pattern = "results_")
species <- gsub("results_","",files)
species <- gsub(".rds","",species)

# Loop over species
library(ggplot2)
library(reshape)
library(truncnorm)

x1 <- seq(-15,20,length.out = 1000)
x2 <- seq(0,520,length.out = 1000)
x3 <- seq(0,35,length.out = 1000)
x4 <- seq(-2,8,length.out = 1000)
x5 <- seq(450,5050,length.out = 1000)
x6 <- seq(0,550,length.out = 1000)
x <- cbind(x1,x2,x3,x4,x5,x6)

burnin <- 10000

# Graphs
i <- 1
for (f in files){
  priors <- readRDS(paste0(dir,files[i]))$priors

  denspriors <- lapply(1:length(names.params), FUN = function(i){
    namedens <- priors[[i]]@distRNG
    substr(namedens,1,1) <- "d"
    do.call(namedens, c(list(x = x[,i]), priors[[i]]@hyperParams))
    })
  denspriors <- do.call(c,denspriors)

  dpriors <- data.frame(x=c(x1,x2,x3,x4,x5,x6),value=denspriors,variable=rep(c("a.cu","b.cu","a.fu","b.fu","mu","s"),each=1000))

  calib <- readRDS(paste0(dir,files[i]))$calib
  head(calib$chain)

  d <- as.data.frame(calib$chain)
  names(d) <- c("a.cu","b.cu","a.fu","b.fu","mu","s")
  row.names(d) <- NULL
  d <- melt(d)
  d$iteration <- rep(1:nrow(calib$chain),6)

  # MCMC traces
  #p <- ggplot(data=d,aes(x=iteration,y=value)) + geom_line() + facet_wrap(~variable, scales = "free_y") + ggtitle(species[i])
  #ggsave(paste0(dir,"mcmc_chain_",species[i],".pdf"),p)

  # posteriors vs priors
  densplot <- ggplot(data=d[-(1:burnin),],aes(x=value)) + geom_histogram(aes(y=..density..),fill="gray60",col="gray60") + facet_wrap(~variable, scales = "free") + ggtitle(species[i])
  densplot <-  densplot + geom_line(data=dpriors,aes(x=x,y=value))
  ggsave(paste0(dir,"posterior_vs_prior_",species[i],".pdf"),densplot)

  # acceptance rates
  ar <- as.data.frame(calib$ar)
  ar <- ar*nrow(calib$ar)
  ar <- lapply(1:6,FUN = function(i){ar[,i]/(1:nrow(ar))})
  ar <- as.data.frame(do.call(cbind,ar))
  names(ar) <- c("a.cu","b.cu","a.fu","b.fu","mu","s")
  ar <- melt(ar)
  ar$iteration <- rep(1:nrow(calib$ar),6)

  #arp <- ggplot(data=ar,aes(x=iteration,y=value)) + geom_line() + facet_wrap(~variable, scales = "free_y") + ggtitle(species[i])
  #ggsave(paste0(dir,"ar_",species[i],".pdf"),arp)

  i <- i+1
}

