# Calibration of the parameters per species
argv <- commandArgs(TRUE)
i <- as.numeric(argv[1])
sizeMC <- as.numeric(argv[2])

devtools::document()
devtools::load_all()

# Import data
files <- list.files(path="data/",pattern = "obs.data")
species <- gsub("obs.data.","",files)
species <- gsub(".rds","",species)

# Parameters shared by all the species
origin.date = "09-01"
var.names = list(date="date",plant="plant",session="session",rep="rep",temp="temp.plant",duration="duration")
temp.params = list(temp.min.cu = -10, temp.max.cu = 15, temp.min.fu = 5, temp.max.fu = 35)
priors = list(a.cu = prior(distRNG="runif", hyperParams=list(min=-5, max=5)),
              b.cu = prior(distRNG="runif", hyperParams=list(min=2, max=200)),
              a.fu = prior(distRNG="runif", hyperParams=list(min=5, max=20)),
              b.fu = prior(distRNG="runif", hyperParams=list(min=1, max=20)),
              mu = prior(distRNG="rnorm", hyperParams=list(mean=1500, sd=1000)),
              s = prior(distRNG="rnorm", hyperParams=list(mean=750, sd=500)))
temp.data <- readRDS("data/temperaturePlants.rds")

# Loop over species

#for (f in files[1]){
obs.data <- readRDS(paste0("data/",files[i]))

calib <- mcmc(data = list(obs.data=obs.data,
                          temp.data=temp.data,
                          var.names=var.names),
              temp.params = temp.params,
              control = list(proposal="AdGl",
                             size = sizeMC))

species <- gsub("obs.data.","",f)
species <- gsub(".rds","",species)

saveRDS(calib,paste0("data/results_",species,".rds"))
