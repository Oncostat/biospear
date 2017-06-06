#################################################################################
#  Survival data tranformation for fitting a Poisson model                      #
#################################################################################
#                                                                               #
#  This function  tranform survival data into a format compatible with          #
#   the glm() function for fitting a Poisson model.                             #
#                                                                               #
#  Its parameters are                                                           #
#   - data          : a data frame with columns:                                #
#                     - id       : the patient identifyier                      #
#                     - time     : the event/censoring time                     #
#                     - status   : the event(1) or censoring(0) indicator       #
#                     - ...      : other factors such like the covariables      #
#                             needed in the regression model                    #
#   - interval.width: the width of the time intervals on wich the risks         #
#                     will be assumed consant                                   #
#   - factors       : a vector of characters, containing the names of the       #
#                     factors to be kept in the transformed data set            #
#   - compress      : a logical, indicating whether the record with the same    #
#                     factor profile should be summarized into one record,      #
#                     i.e. whether the data should be expressed in a short form #
#                                                                               #
#################################################################################
#                                                                               #
# Author: Federico Rotolo <federico.rotolo@gustaveroussy.fr>                    #
# Original code by Stephanie Kovalchik                                          #
# http://r.789695.n4.nabble.com/exponential-proportional-hazard-model-td805536.html #
#                                                                               #
#   Date:                 June 22, 2015                                         #
#   Last modification on: July 08, 2015                                         #
#################################################################################
poissonize <- function(data, interval.width = 365.25/4,
                       factors = NULL, compress = TRUE) {
  # THE FOLLOWING FUNCTION COMPUTES ALL THE TIME BREAKS GIVEN THE EXIT TIME
  person.breaks <- Vectorize(function(stopT)
    unique(c(seq(0, stopT, by = interval.width), stopT)))
  
  # NEXT WE GET A LIST OF EACH SUBJECT'S TIME PERIODS
  the.breaks <- person.breaks(data$time)
  
  # NOW WE OBTAIN THE EXPANDED PERIOD OF OBSERVATION
  startT <- lapply(the.breaks, function(x) x[-length(x)])  # LEFT TIME POINT
  stopT <-  lapply(the.breaks, function(x) x[-1])    # RIGHT TIME POINTS
  
  # THE FOLLOWING ARE NEEDED TO COMPLETE THE LONG VERSION OF THE DATA SET
  count.per.id <- sapply(startT, length)
  index <- tapply(data$id, data$id, length)
  index <- cumsum(index) # INDEX OF LAST OBSERVATION FOR EACH PATIENT
  event <- rep(0, sum(count.per.id))
  event[cumsum(count.per.id)] <- data$status[index]
  
  # BRING ALL OF THIS TOGETHER TO CREATE THE EXPANDED DATASET
  expData <- cbind(
    data.frame(id = rep(data$id[index], count.per.id),
               startT = unlist(startT),
               stopT = unlist(stopT),
               event = event),
    data[sapply(rep(data$id[index], count.per.id), 
                function(x) which(data$id==x)), factors, drop=FALSE])
  
  
  
  # CREATE TIME VARIABLE WHICH INDICATES THE PERIOD OF OBSERVATION
  # THIS WILL BE THE OFFSET IN THE POISSON MODEL FIT
  expData$time <- expData$stopT - expData$startT # LENGTH OF OBSERVATION
  
  # NEXT WE CREATE A FACTOR FOR EACH INTERVAL THAT WILL ALLOW US
  # TO HAVE A DIFFERENT RATE FOR EACH PERIOD
  expData$interval <- factor(expData$start)
  expData <- subset(expData, select=-c(startT, stopT))
  expData$time <- expData$time
  
  if (compress) {
    m <- aggregate(x = expData$event,
                   by = expData[, c('interval', factors), drop=FALSE], 
                   sum)
    colnames(m)[ncol(m)] <- 'm'
    Rt <- aggregate(x = expData$time,
                    by = expData[, c('interval', factors), drop=FALSE], 
                    sum)
    colnames(Rt)[ncol(Rt)] <- 'Rt'
    N <- aggregate(x = expData$event,
                   by = expData[, c('interval', factors), drop=FALSE], 
                   function(x) sum(!is.na(x)))
    colnames(N)[ncol(N)] <- 'N'
    expData <- base::merge(m, base::merge(Rt, N))
  }
  
  attr(expData, 'interval.width') <- interval.width
  return(expData)
}
################################################################################
