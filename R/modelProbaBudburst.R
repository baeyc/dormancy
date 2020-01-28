#' Model for probability of budburst
#'
#' Compute the probability of budburst as a function of accumulated CU and Fu
#'
#' @name modelProbaBB
#' @aliases modelProbaBB
#'
#' @param temp.data a data frame with the extracted temperatures experienced by the plant (as a result of the \code{\link{extractTemp}} function)
#' @param var.names a list with the names of the variables containing the date, the plant ID, the session, the repetition and the temperature. Default is
#' \code{var.names = list(date="date",plant="plant",session="session",rep="rep",temp="temp.plant")}
#' @param temp.param a list of min and max temperatures for the CU and FU functions, containing \code{temp.min.cu, temp.min.fu,temp.max.cu, temp.max.fu}
#' @param cufu.params a list of parameters for the CU and FU curves. It must contain the following parameters: \code{a.cu, b.cu, a.fu, b.fu, mu, s}
#' @param origin.date the origin date when plants start accumulating temperatures, in the format "mm-dd" (e.g. "09-01" for 1st of September)
#' @return the probability of budburst per plant, repetition and session
#'
#' @export modelProbaBB
modelProbaBB <- function(temp.data,
                         var.names = list(date="date",plant="plant",session="session",rep="rep",temp="temp.plant"),
                         temp.params,
                         cufu.params,origin.date="09-01"){
  # create a copy of temp.data
  data <- temp.data

  # compute CU and FU
  data$cu <- chillingUnits(data,
                           var.names = list(temp=var.names$temp,date=var.names$date),
                           temp.min = temp.params$temp.min.cu, temp.max= temp.params$temp.max.cu,
                           mu = cufu.params$a.cu, s = cufu.params$b.cu)
  data$fu <- forcingUnits(data,
                          var.names = list(temp=var.names$temp,date=var.names$date),
                          temp.min = temp.params$temp.min.fu, temp.max= temp.params$temp.max.fu,
                          a = cufu.params$a.fu, b = cufu.params$b.fu)

  # compute accumulation from origin date
  # add origin
  year <- lubridate::year(data[,var.names$date])
  dateOrigin <- paste0(year,"-",origin.date)
  doyOrigin <- lubridate::yday(lubridate::as_date(dateOrigin))
  data$before.origin <- ifelse(lubridate::yday(data$date) < doyOrigin, TRUE, FALSE)
  # do not account for temperatures accumulated between 01/01 and origin.date
  data$fu[data$before.origin] <- 0
  data$cu[data$before.origin] <- 0

  data <- data %>% dplyr::group_by_at(c(var.names$session,var.names$plant,var.names$rep)) %>% dplyr::mutate(cu.cum = cumsum(cu), fu.cum = cumsum(fu))

  dataPerPlantSessionRep <- data %>%
    dplyr::group_by_at(c(var.names$session,var.names$plant,var.names$rep)) %>%
    dplyr::mutate(cu = max(cu.cum), fu = max(fu.cum)) %>%
    dplyr::select(var.names$session,var.names$plant,var.names$rep,cu,fu) %>%
    dplyr::distinct()

  dataPerPlantSessionRep$probaBB <- 1/(1+exp(-(dataPerPlantSessionRep$cu + dataPerPlantSessionRep$fu - cufu.params$mu)/cufu.params$s))

  return(dataPerPlantSessionRep)
}
