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
temp.params = list(temp.min.cu = -10, temp.max.cu = 15, temp.min.fu = 0, temp.max.fu = 35)
priors = list(a.cu = prior(distRNG="rnorm", hyperParams=list(mean=7, sd=3)),
              b.cu = prior(distRNG="rtruncnorm", hyperParams=list(a=2, mean=100, sd=40)),
              a.fu = prior(distRNG="rnorm", hyperParams=list(mean=5, sd=4)),
              b.fu = prior(distRNG="rnorm", hyperParams=list(mean=2.5, sd=1)),
              mu = prior(distRNG="rnorm", hyperParams=list(mean=500, sd=150)),
              s = prior(distRNG="rtruncnorm", hyperParams=list(a=0,mean=50, sd=20)))
temp.data <- readRDS("data/temperaturePlants.rds")

# # Aggregate temperatures: one for the day and one for the night
# temp.data$day <- as.Date(temp.data$date) # remove hours
# temp.aggregated <- temp.data %>% dplyr::group_by(day) %>% dplyr::mutate(temp.day = min(temp.plant), temp.night = max(temp.plant)) # compute min and max per day
# temp.data <- temp.aggregated %>% dplyr::select(-date,-temp.out,-temp.in,-temp.plant) %>% dplyr::distinct() # drop variables at the 3-hour level
#
# # create a unique column from cols "temp.day" and "temp.night"
# temp.data <- tidyr::pivot_longer(temp.data, cols = tidyr::starts_with("temp"), values_to = "temp.plant") %>% dplyr::select(-name)

# Add species in temp.data
temp.data$species <- gsub('[0-9]+', '', temp.data$plant)

# Loop over species

#for (f in files[1]){
obs.data <- readRDS(paste0("data/",files[i]))

calib <- mcmc(data = list(obs.data=obs.data,
                          temp.data=temp.data,
                          var.names=var.names),
              temp.params = temp.params,
	      priors = priors,
              control = list(proposal="AdGl",
                             size = sizeMC))

species <- gsub("obs.data.","",files[i])
species <- gsub(".rds","",species)

saveRDS(list(calib=calib,priors=priors),paste0("data/results_",species,".rds"))


