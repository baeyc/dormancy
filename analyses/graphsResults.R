# Analyze MCMC results
setwd("data/results/AdaptiveGlobal/")

# Import data
files <- list.files(pattern = "results_")
species <- gsub("results_","",files)
species <- gsub(".rds","",species)

# Loop over species
library(ggplot2)
library(reshape)


x1 <- seq(-5,5,length.out = 1000)
x2 <- seq(2,200,length.out = 1000)
x3 <- seq(5,20,length.out = 1000)
x4 <- seq(1,20,length.out = 1000)
x5 <- seq(-2000,4000,length.out = 1000)
x6 <- seq(-1000,3000,length.out = 1000)
priors <- c(dunif(x1,-5,5),
            dunif(x2,2,200),
            dunif(x3,5,20),
            dunif(x4,1,20),
            dnorm(x5,1500,1000),
            dnorm(x6,750,500))

dpriors <- data.frame(x=c(x1,x2,x3,x4,x5,x6),value=priors,variable=rep(c("a.cu","b.cu","a.fu","b.fu","mu","s"),each=1000))


# Graphs
i <- 1
for (f in files){
  calib <- readRDS(f)
  head(calib$chain)

  d <- as.data.frame(calib$chain)
  names(d) <- c("a.cu","b.cu","a.fu","b.fu","mu","s")
  row.names(d) <- NULL
  d <- melt(d)
  d$iteration <- rep(1:nrow(calib$chain),6)

  # MCMC traces
  p <- ggplot(data=d,aes(x=iteration,y=value)) + geom_line() + facet_wrap(~variable, scales = "free_y") + ggtitle(species[i])
  ggsave(paste0("mcmc_chain_",species[i],".pdf"),p)

  # posteriors vs priors
  densplot <- ggplot(data=d[-(1:burnin),],aes(x=value)) + geom_histogram(aes(y=..density..),fill="gray60",col="gray60") + facet_wrap(~variable, scales = "free") + ggtitle(species[i])
  densplot <-  densplot + geom_line(data=dpriors,aes(x=x,y=value))
  ggsave(paste0("posterior_vs_prior_",species[i],".pdf"),densplot)

  # acceptance rates
  ar <- as.data.frame(calib$ar)
  ar <- ar*nrow(calib$ar)
  ar <- lapply(1:6,FUN = function(i){ar[,i]/(1:nrow(ar))})
  ar <- as.data.frame(do.call(cbind,ar))
  names(ar) <- c("a.cu","b.cu","a.fu","b.fu","mu","s")
  ar <- melt(ar)
  ar$iteration <- rep(1:nrow(calib$ar),6)

  arp <- ggplot(data=ar,aes(x=iteration,y=value)) + geom_line() + facet_wrap(~variable, scales = "free_y") + ggtitle(species[i])
  ggsave(paste0("ar_",species[i],".pdf"),arp)

  i <- i+1
}
