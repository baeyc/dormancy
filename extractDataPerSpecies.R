# Example
devtools::document()
devtools::load_all()

# Import weather data
metlille <- read.table("data/meteo.Lesquin.txt", sep="\t", header=T) # in Lesquin
metlille$mydate <- apply(metlille, MARGIN=1, FUN=function(li) paste(li["year"], li["month"], li["day"], li["hour"], sep=" "))
metlille$posixct <- as.POSIXct(metlille$mydate, format="%Y %m %d %H", tz="GMT")
metlille$temp <- metlille$t - 273.15 # from Kelvin to Celsius

temp.outside <- metlille
var.names.out <- list(temp="temp",date="posixct")

metser <- read.table("data/meteo.serre", sep="\t", header=T) # in the greenhouse
metser$Date <- as.POSIXct(metser$Date, format="%d.%m.%Y %H:%M", tz="GMT")
metsersub <- metser[seq(1,dim(metser)[1],15),] # every 3 hours
metsersub$Date <- metsersub$Date - lubridate::minutes(6) # go back to round hours

temp.inside <- metsersub
var.names.in <- list(temp="Tint",date="Date")

# Import plant data
bb <- read.table("data/mydata.bb.txt", sep="\t", header=T)
bb$date_collecte <- as.POSIXct(bb$date_collecte, format="%d/%m/%Y")
bb$date_debourrement <- as.POSIXct(bb$date_debourrement, format="%d/%m/%Y")

data.deb <- bb
var.names.deb <- list(session="session",plant="Individu",harv="date_collecte",budburst="date_debourrement")

# Extract temperatures experienced by each plant
temp.plants <- extractTemp(temp.outside,temp.inside,var.names.out,var.names.in,data.deb,var.names.deb)

# add origin
origin.date <- "09-01"
year <- lubridate::year(temp.plants$date)
dateOrigin <- paste0(year,"-",origin.date)
doyOrigin <- lubridate::yday(lubridate::as_date(dateOrigin))
temp.plants$before.origin <- ifelse(lubridate::yday(temp.plants$date) < doyOrigin, TRUE, FALSE)

# get elapsing time between two successive temperature recording
duration <- as.double(abs(difftime(temp.plants$date,dplyr::lead(temp.plants$date),units="hours"))) # time in hours elapsing between two successive recordings
temp.plants$duration <- duration

saveRDS(temp.plants,"data/temperaturePlants.rds")


# Compute the probability of budburst for each plant, as a function of accumulated CU and FU
probaBB <- modelProbaBB(temp.data = temp.plants,
                        origin.date = "09-01",
                        var.names = list(date="date",plant="plant",session="session",rep="rep",temp="temp.plant"),
                        temp.params = list(temp.min.cu = -10, temp.max.cu = 15, temp.min.fu = 5, temp.max.fu = 35),
                        cufu.params = list(a.cu = 5, b.cu = 5, a.fu = 15, b.fu = 4, mu = 1500, s = 750))

# get rep ID for each harvesting date
rep.session <- temp.plants %>%
  dplyr::select(session,plant,rep,harv.date) %>%
  dplyr::distinct()

rep.session <- dplyr::inner_join(rep.session, data.deb, by = c("session"="session","plant"="Individu", "harv.date"="date_collecte")) %>%
  dplyr::select(session,plant,Espece,rep,harv.date,date_debourrement)
head(rep.session)
rep.session$budburst <- 1*!is.na(rep.session$date_debourrement)
rep.session$date_debourrement <- NULL


# Divide the data per plant species (each species will be treated separately)
species <- unique(rep.session$Espece)
data.per.species <- lapply(species,FUN = function(s){rep.session[rep.session$Espece==s,]})
names(data.per.species) <- species

lapply(species,FUN=function(s){saveRDS(data.per.species[[s]],paste0("data/obs.data.",s,".rds"))})

