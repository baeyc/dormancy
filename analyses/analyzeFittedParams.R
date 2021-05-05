# Analyze MCMC results
devtools::document()
devtools::load_all()

# Import data
dir <- "data/results/AdaptiveGlobal/Priors unif et norm/"
dirArticle <- "~/Documents/Articles/En cours/Dormancy/figures/"
files <- list.files(path = dir, pattern = "results_")
species <- gsub("results_","",files)
species <- gsub(".rds","",species)
species_long <- c("Acer pseudoplatanus", "Aesculus hippocastanum", "Alnus glutinosa",
                  "Betula pendula", "Carpinus betulus", "Castanea sativa", "Corylus avellana",
                  "Fagus sylvativa", "Fraexinus excelsior", "Larix decidua", "Prunus avium", "Quercus robur")


# Loop over species
library(ggplot2)
library(plotROC)
library(reshape)

burnin <- 10000

# Posterior mean
postmean <- lapply(1:length(files), FUN = function(i){
  calib <- readRDS(paste0(dir,files[i]))$calib
  chain_after_burnin <- calib$chain[-(1:burnin),]
  d <- as.data.frame(t(apply(chain_after_burnin,2,mean)))
  names(d) <- c("a.cu","b.cu","a.fu","b.fu","mu","s")
  d$species <- species[i]
  d
})
postmean <- do.call(rbind,postmean)

# Maximum a posteriori
map <- lapply(1:length(files), FUN = function(i){
  calib <- readRDS(paste0(dir,files[i]))$calib
  chain_after_burnin <- calib$chain[-(1:burnin),]
  d <- apply(chain_after_burnin,2,density)
  map <- sapply(1:length(d),FUN=function(j){
    d[[j]]$x[which.max(d[[j]]$y)]
  })
  map <- as.data.frame(t(map))
  names(map) <- c("a.cu","b.cu","a.fu","b.fu","mu","s")
  map$species <- species[i]
  map
})
map <- do.call(rbind,map)

# CU and FU curves
a1 <- (map$a.cu+10)/40 * (map$b.cu-2) + 1
b1 <- (map$b.cu-2)*(1-(map$a.cu+10)/40) + 1
a2 <- (postmean$a.cu+10)/40 * (postmean$b.cu-2) + 1
b2 <- (postmean$b.cu-2)*(1-(postmean$a.cu+10)/40) + 1

x <- seq(-10,30,0.05)
xnorm <- (x-(-10))/40

y <- numeric()
y2 <- numeric()
for (i in 1:6)
{
  mode <- (a1[i]-1)/(a1[i]+b1[i]-2)
  max <- dbeta(mode,a1[i],b1[i])
  y <- c(y,dbeta(xnorm,a1[i],b1[i])/max)
  mode <- (a2[i]-1)/(a2[i]+b2[i]-2)
  max <- dbeta(mode,a2[i],b2[i])
  y2 <- c(y2,dbeta(xnorm,a2[i],b2[i])/max)
}

d <- data.frame(x=rep(x,12),y=y,a=rep(a1,each=length(x)),b=rep(b1,each=length(x)),species=rep(species,each=length(x)))
d2 <- data.frame(x=rep(x,12),y=y2,a=rep(a2,each=length(x)),b=rep(b2,each=length(x)),species=rep(species,each=length(x)))

p <- ggplot(data=d,aes(x=x,y=y,color=species)) + geom_line(lwd=1.5) + scale_color_discrete(name = "Species") + xlab("Temperature") + ylab("CU")
p <- p + xlim(c(0,20)) + facet_wrap(~species,ncol=4)
ggsave(paste0(dir,"cu_curves_map.pdf"),p)

p2 <- ggplot(data=d2,aes(x=x,y=y,color=species)) + geom_line(lwd=1.5) + scale_color_discrete(name = "Species") + xlab("Temperature") + ylab("CU")
p2 <- p2 + xlim(c(0,20)) + facet_wrap(~species,ncol=4)

dall <- rbind(d,d2)
dall$type <- rep(c("MAP","Posterior mean"),each=nrow(d))
dall$species <- factor(dall$species, levels = species, labels = species_long)

p3 <- ggplot(data=dall,aes(x=x,y=y,color=species,linetype=type)) + geom_line() + scale_color_discrete(name = "Species") + xlab("Temperature") + ylab("CU")
p3 <- p3 + xlim(c(0,20)) + facet_wrap(~species,ncol=4) + scale_linetype_discrete("Estimation") + guides(color=FALSE) + theme(strip.text = element_text(face = "italic"))
p3
ggsave(paste0(dirArticle,"cu_curves.pdf"),p3)


## FU
a1 <- map$a.fu
b1 <- map$b.fu
a2 <- postmean$a.fu
b2 <- postmean$b.fu

x <- seq(0,35,0.05)

y <- numeric()
y2 <- numeric()
for (i in 1:length(a1))
{
  y <- c(y,1/(1+exp(-(x-a1[i])/b1[i])))
  y2 <- c(y2,1/(1+exp(-(x-a2[i])/b2[i])))
}

d <- data.frame(x=rep(x,length(a1)),y=y,c=rep(a1,each=length(x)),d=rep(b1,each=length(x)),species=rep(species,each=length(x)))
d2 <- data.frame(x=rep(x,length(a1)),y=y2,c=rep(a2,each=length(x)),d=rep(b2,each=length(x)),species=rep(species,each=length(x)))
dall <- rbind(d,d2)
dall$type <- rep(c("MAP","Posterior mean"),each=nrow(d))
dall$species <- factor(dall$species, levels = species, labels = species_long)

p <- ggplot(data=d,aes(x=x,y=y,color=species)) + geom_line() + scale_color_discrete(name = "Species") + xlab("Temperature") + ylab("FU") + geom_vline(xintercept = 5,linetype=2,col="darkgrey") +
  geom_vline(xintercept = 35,linetype=2,col="darkgrey")  + facet_wrap(~species,ncol=4)
p

pp <- ggplot(data=dall,aes(x=x,y=y,color=species,linetype=type)) + geom_line() + guides(color=FALSE) + xlab("Temperature") + ylab("FU") + geom_vline(xintercept = 5,linetype=2,col="darkgrey") +
  geom_vline(xintercept = 35,linetype=2,col="darkgrey")  + facet_wrap(~species,ncol=4) + scale_linetype_discrete(name="Estimation")
pp <- pp + theme(strip.text = element_text(face = "italic"))
ggsave(paste0(dirArticle,"fu_curves.pdf"),pp)
