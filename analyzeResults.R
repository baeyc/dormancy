# Analyze MCMC results
devtools::document()
devtools::load_all()

# Import data
files <- list.files(pattern = "results_")
species <- gsub("results_","",files)
species <- gsub(".rds","",species)

# Loop over species
library(ggplot2)
library(plotROC)
library(reshape)

burnin <- 10000

# Posterior mean
postmean <- lapply(1:length(files), FUN = function(i){
  calib <- readRDS(paste0("data/results/AdaptiveGlobal/",files[i]))
  chain_after_burnin <- calib$chain[-(1:burnin),]
  d <- as.data.frame(t(apply(chain_after_burnin,2,mean)))
  names(d) <- c("a.cu","b.cu","a.fu","b.fu","mu","s")
  d$species <- species[i]
  d
})
postmean <- do.call(rbind,postmean)

# Maximum a posteriori
map <- lapply(1:length(files), FUN = function(i){
  calib <- readRDS(paste0("data/results/AdaptiveGlobal/",files[i]))
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


# Compare predictions using posterior means and MAP with observations
temp.data <- readRDS("data/temperaturePlants.rds")
origin.date = "09-01"
var.names = list(date="date",plant="plant",session="session",rep="rep",temp="temp.plant",duration="duration")
temp.params = list(temp.min.cu = -10, temp.max.cu = 15, temp.min.fu = 5, temp.max.fu = 35)
names.params <- c("a.cu","b.cu","a.fu","b.fu","mu","s")

obs.files <- list.files(path="data/",pattern = "obs.data")

pred <- lapply(1:length(files), FUN = function(i){
  params1 <- postmean[i,]
  params2 <- map[i,]
  obs.data <- readRDS(paste0("data/",obs.files[i]))

  predPostMean <- modelProbaBB(temp.data = temp.data,
                      var.names = var.names,
                      temp.params = temp.params,
                      cufu.params =stats::setNames(as.list(params1),names.params))
  names(predPostMean) <- c("session","plant","rep","cu.postmean","fu.postmean","probaBB.postmean")
  predMAP <- modelProbaBB(temp.data = temp.data,
                               var.names = var.names,
                               temp.params = temp.params,
                               cufu.params =stats::setNames(as.list(params2),names.params))
  names(predMAP) <- c("session","plant","rep","cu.map","fu.map","probaBB.map")

  pred <- dplyr::inner_join(predMAP,predPostMean,by = c("session", "plant", "rep"))
  pred <- dplyr::inner_join(obs.data,pred,by = c("session", "plant", "rep"))

  dpred <- melt(pred,id.vars = "budburst", measure.vars = c("probaBB.postmean","probaBB.map"))
  #ggplot(data=dpred,aes(x=as.factor(budburst),y=value,fill=variable)) + geom_boxplot() + xlab("Observed budburst") + ylab("Predicted probability of budburst") +
  #  scale_fill_discrete(name="Estimate",labels=c("posterior mean","maximum a posteriori"))
  #ggsave(paste0("data/results/AdaptiveGlobal/predictions_",species[i],".pdf"),height=5,width=6)

  roc <- ggplot(dpred,aes(m=value,d=budburst,color=variable)) + geom_roc(labels=F,pointsize=0)
  auc <- calc_auc(roc)
  col <- unique(ggplot_build(roc)$data[[1]]["colour"])
  roc + style_roc(theme = theme_gray,xlab = "1 - Specificity") + scale_color_discrete(name="Estimate",labels=c("posterior mean","maximum a posteriori")) +
    annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(roc)$AUC[1], 2)), color=col$colour[1]) +
    annotate("text", x = .75, y = .2, label = paste("AUC =", round(calc_auc(roc)$AUC[2], 2)), color=col$colour[2])
  ggsave(paste0("data/results/AdaptiveGlobal/roc_",species[i],".pdf"),height=5,width=6)
})



