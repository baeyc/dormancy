# Analyze MCMC results
devtools::document()
devtools::load_all()

# Import data
dir <- "data/results/AdaptiveGlobal/Priors unif et norm/"
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

ci95 <- lapply(1:length(files), FUN = function(i){
  calib <- readRDS(paste0(dir,files[i]))$calib
  chain_after_burnin <- calib$chain[-(1:burnin),]
  d <- as.data.frame(t(apply(chain_after_burnin,2,quantile,c(0.025,0.975))))
  d$param <- c("a.cu","b.cu","a.fu","b.fu","mu","s")
  d$species <- species[i]
  d
})
ci95 <- do.call(rbind,ci95)

# Table with MAP and 95% CI
map_t <- melt(map)
names(map_t) <- c("species","param","MAP")
map_ci95 <- merge(map_t,ci95,by=c("species","param"))

ci952 <- melt(ci95)
ci952$var <- paste0(ci952$param,ci952$variable)
ci952$param <- ci952$variable <- NULL
ci952 <- dcast(ci952, species ~ var)
ci95merged <- ci952$species
for (i in 1:6){
  ci95merged <- cbind(ci95merged,paste0("[",format(ci952[,2*i],digits=3)," ; ",format(ci952[,2*i+1],digits=3),"]"))
}
ci95merged <- as.data.frame(ci95merged)
names(ci95merged) <- c("species",paste0(c("a.cu","b.cu","a.fu","b.fu","mu","s"),"_95CI"))

map_ci95_2 <- merge(map,ci95merged)

library(xtable)
print(xtable(map_ci95_2[,c(1,2,8,3,9,4,10,5,11)]), include.rownames=F)
print(xtable(map_ci95_2[,c(1,6,12,7,13)]), include.rownames=F)

# Compare predictions using posterior means and MAP with observations
temp.data <- readRDS("data/temperaturePlants.rds")
origin.date = "09-01"
var.names = list(date="date",plant="plant",session="session",rep="rep",temp="temp.plant",duration="duration")
temp.params = list(temp.min.cu = -10, temp.max.cu = 15, temp.min.fu = 5, temp.max.fu = 35)
names.params <- c("a.cu","b.cu","a.fu","b.fu","mu","s")

obs.files <- list.files(path="data/",pattern = "obs.data")

dpredall <- lapply(1:length(files), FUN = function(i){
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
  dpred$species <- species[i]
  dpred$n <- nrow(pred)
  # pred_plot <- ggplot(data=dpred[dpred$variable=="probaBB.map",],aes(x=as.factor(budburst),y=value)) +
  #   geom_boxplot(width=0.25,fill="grey") + xlab("Observed budburst") + ylab("Predicted probability of budburst") +
  #   #scale_fill_discrete(name="Estimate",labels=c("posterior mean","maximum a posteriori")) +
  #   guides(fill = FALSE) +
  #   ggtitle(paste0(species_long[i],", n=",nrow(pred))) + ylim(c(0,1))
  # #ggsave(paste0(dir,"predictions_",species[i],".pdf"),pred_plot,height=5,width=6)
  #
  # roc <- ggplot(dpred[dpred$variable=="probaBB.map",],aes(m=value,d=budburst)) + geom_roc(labels=F,pointsize=0) + ggtitle(paste0(species_long[i],", n=",nrow(pred)))
  # auc <- calc_auc(roc)
  # col <- unique(ggplot_build(roc)$data[[1]]["colour"])
  # roc <- roc + style_roc(theme = theme_gray,xlab = "1 - Specificity") +
  #   #scale_color_discrete(name="Estimate",labels=c("posterior mean","maximum a posteriori")) +
  #   guides(fill = FALSE) +
  #   annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(roc)$AUC[1], 2)), color=col$colour[1]) +
  #   annotate("text", x = .75, y = .2, label = paste("AUC =", round(calc_auc(roc)$AUC[2], 2)), color=col$colour[2])
  # #ggsave(paste0(dir,"roc_",species[i],".pdf"),roc, height=5,width=6)

  #return(list(pred=pred_plot,roc=roc))
  return(dpred)
})

dpredall <- do.call(rbind,dpredall)
dpredall$species <- factor(dpredall$species, levels = species, labels = species_long)
dpredall$title <- paste0(dpredall$species,", n=",dpredall$n)

pred_plot <- ggplot(data=dpredall[dpredall$variable=="probaBB.map",],aes(x=as.factor(budburst),y=value)) +
     geom_boxplot(width=0.25,fill="grey") + xlab("Observed budburst") + ylab("Predicted probability of budburst") +
     #scale_fill_discrete(name="Estimate",labels=c("posterior mean","maximum a posteriori")) +
     guides(fill = FALSE) + facet_wrap(~title,nc=4,nr=3) +
     ylim(c(0,1)) + theme(strip.text = element_text(face = "italic"))
pred_plot
ggsave("~/Documents/Articles/En cours/Dormancy/figures/pred_plot.pdf",pred_plot,height = 7, width = 8.5)

rocf <- ggplot(dpredall[dpredall$variable=="probaBB.map",],aes(m=value,d=budburst)) + geom_roc(labels=F,pointsize=0,size=0.5) +
  facet_wrap(~title,nc=4,nr=3)
auc <- calc_auc(rocf)
auc$title <- unique(dpredall$title)
#col <- unique(ggplot_build(roc)$data[[1]]["colour"])
roc <- rocf + style_roc(theme = theme_gray,xlab = "1 - Specificity") +
 #scale_color_discrete(name="Estimate",labels=c("posterior mean","maximum a posteriori")) +
 guides(fill = FALSE) + #facet_wrap(~title,nc=3,nr=4) +
 geom_text(data = auc, aes(label=paste0("AUC=",round(AUC,2))),
             x = Inf, y = -Inf, hjust=1, vjust=0,
             inherit.aes = FALSE) + theme(strip.text = element_text(face = "italic"))
 #annotate("text", x = 0.75, y=0.25, label = auc$AUC, group=auc$AUC)
 #annotate("text", x = .75, y = .25, label = paste("AUC =", auc$AUC))#, color=col$colour[1])
 #annotate("text", x = .75, y = .2, label = paste("AUC =", round(calc_auc(roc)$AUC[2], 2)), color=col$colour[2])
roc
ggsave("~/Documents/Articles/En cours/Dormancy/figures/roc_all.pdf",roc,height = 7, width = 8.5)


# Plot estimated CU+FU as a function of time
lapply(1:length(files), FUN = function(i){
  params1 <- postmean[i,]
  params2 <- map[i,]
  obs.data <- readRDS(paste0("data/",obs.files[i]))

  data <- dplyr::right_join(temp.data,obs.data)

  # compute CU and FU
  data$cu.mean <- chillingUnits(data,
                                var.names = var.names,
                                temp.min = temp.params$temp.min.cu, temp.max= temp.params$temp.max.cu,
                                mu = params1$a.cu, s = params1$b.cu)
  data$fu.mean <- forcingUnits(data,
                               var.names = var.names,
                               temp.min = temp.params$temp.min.fu, temp.max= temp.params$temp.max.fu,
                               a = params1$a.fu, b = params1$b.fu)

  # do not account for temperatures accumulated between 01/01 and origin.date
  data$fu.mean[data$before.origin] <- 0
  data$cu.mean[data$before.origin] <- 0

  data <- data %>% dplyr::group_by_at(c(var.names$session,var.names$plant,var.names$rep)) %>% dplyr::mutate(cu.cum = cumsum(cu.mean), fu.cum = cumsum(fu.mean))
  dataCU <- data %>% dplyr::select(!fu.cum)
  dataCU$units <- dataCU$cu.cum
  dataCU$type <- "Chilling"
  dataFU <- data %>% dplyr::select(!cu.cum)
  dataFU$units <- dataFU$fu.cum
  dataFU$type <- "Forcing"
  dataCU$cu.cum <- NULL
  dataFU$fu.cum <- NULL

  dataCUFU <- rbind(dataCU,dataFU)

  dataBB <- data %>% dplyr::select(session,rep,plant,harv.date,budburst) %>% dplyr::distinct()

  p1 <- ggplot(data,aes(x=date,y=cu.cum,col=as.factor(rep)),group=as.factor(rep)) + geom_line() +
    facet_wrap(~session,scales = "free") + scale_color_discrete(name="Branch")

  p12 <- ggplot(dataCUFU,aes(x=date,y=units,col=as.factor(rep),linetype=type),group=as.factor(rep)) + geom_line() +
    facet_wrap(~session,scales = "free") + scale_color_discrete(name="Branch") + scale_linetype_discrete(name="Units")

  p2 <- ggplot(data,aes(x=date,y=fu.cum,col=as.factor(rep)),group=as.factor(rep)) + geom_line() +
    facet_wrap(~session,scales = "free") + scale_color_discrete(name="Branch")
  ggsave(paste0(dir,"/cu.cum_",species[i],".pdf"),p1,height = 6, width = 12)
  ggsave(paste0(dir,"fu.cum_",species[i],".pdf"),p2,height = 6, width = 12)
  ggsave(paste0(dir,"cufu.cum_",species[i],".pdf"),p12,height = 6, width = 12)

  data$probaBB <- 1/(1+exp(-(data$cu.cum + data$fu.cum - params1$mu)/params1$s))

  dataBB <- data %>% dplyr::select(session,harv.date,plant,rep,budburst) %>% dplyr::distinct()

  p3 <- ggplot(data,aes(x=date,y=probaBB,col=as.factor(rep)),group=as.factor(rep)) + geom_line() + facet_grid(session~rep,scales = "free", labeller=label_both) + ylim(c(0,1))
  p3 <- p3 + geom_point(data=dataBB,aes(x=harv.date,y=budburst,pch='.')) + scale_color_discrete(name="Branch") + scale_shape_discrete(name="Budburst") + theme(axis.text.x=element_text(angle=45, hjust=1))

  ggsave(paste0(dir,"probaBB_",species[i],".pdf"),p3,height = 6, width = 12)
})
