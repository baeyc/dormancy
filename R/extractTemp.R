#' Temperatures per plant
#'
#' Extract the temperatures experienced by each plant, taking into account both periods before and after harvest
#'
#' @name extractTemp
#' @aliases extractTemp
#'
#' @param temp.outside a data frame with the dates and temperatures outside the greenhouse
#' @param temp.inside a data frame with the dates and temperatures inside the greenhouse
#' @param var.names.out the names of the temperature and date variables, in the format \code{list(temp="temp.name",date="date.name")} for \code{temp.outside}
#' @param var.names.in the names of the temperature and date variables, in the format \code{list(temp="temp.name",date="date.name")} for \code{temp.inside}
#' @param data.deb a data frame with the information about budburst
#' @param var.names.deb the names of the variables representing the session, the plant ID, the harvesting date and the budburst status in \code{data.deb},
#' in the format \code{list(session="sess.name",plant="plant.name",harv="harv.name",budburst="bud.name")}
#' @return a data frame with the accumulated temperatures outside and inside the greenhouse by the plant from 1st of january to 31st of december
#'
#' @export extractTemp

extractTemp <- function(temp.outside,temp.inside,var.names.out,var.names.in,
                       data.deb,var.names.deb){


  # data frame final : session, individu, espèce, jour, outside(0/1), deb
  df <- numeric()
  plants <- unique(data.deb[,var.names.deb$plant]) # get unique IDs of plant

  perPlant <- lapply(plants, FUN = function(i){
    sessions <- unique(data.deb[data.deb[,var.names.deb$plant]==i,var.names.deb$session]) # get the sessions where plant i was sampled

    perSession <- lapply(sessions, FUN = function(k){
      # get the years of session k and create a subset of the original data corresponding to plant i, session k and years of that session
      yearsInSession <- unique(lubridate::year(data.deb[data.deb[,var.names.deb$plant]==i & data.deb[,var.names.deb$session]==k,var.names.deb$harv]))
      temp.out.session <- temp.outside[lubridate::year(temp.outside[,var.names.out$date])%in%yearsInSession,]
      temp.in.session <- temp.inside[lubridate::year(temp.inside[,var.names.in$date])%in%yearsInSession,]
      data.deb.session <- data.deb[data.deb[,var.names.deb$session]==k & data.deb[,var.names.deb$plant]==i,]

      # get the dates at which data were collected for plant i at sessions k
      time.meas <- unique(data.deb.session[,var.names.deb$harv])

      perRep <- lapply(1:length(time.meas), FUN= function(j){
        # create a data frame with the sequence of dates of the session and the temperatures experienced by plant (from temp.outside before harvest
        # and from temp.indise after harvest)
        time.seq <- seq(ISOdate(min(yearsInSession),1,1,0,0,0,tz="GMT"), ISOdate(max(yearsInSession),04,30,tz="GMT"), "3 hours")
        d <- data.frame(session=k,plant=i,date=time.seq,rep=j,harv.date=time.meas[j],outside=(time.seq<time.meas[j]))

        # merge with temperature data
        d <- merge(d,temp.out.session[,c(var.names.out$temp,var.names.out$date)],by.x=c("date"),by.y=var.names.out$date)
        names(d)[names(d)==var.names.out$temp] <- "temp.out"
        d <- merge(d,temp.in.session[,c(var.names.in$temp,var.names.in$date)],by.x=c("date"),by.y=var.names.in$date)
        names(d)[names(d)==var.names.in$temp] <- "temp.in"
        d$temp.plant <- ifelse(d$outside,d$temp.out,d$temp.in)

        df <- rbind(df,d)
        return(df)
      })
      perRep <- do.call(rbind,perRep)
      return(perRep)
    })
    perSession <- do.call(rbind,perSession)
    return(perSession)
  })
  df <- do.call(rbind,perPlant)

  return(df)
}
