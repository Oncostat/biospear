################################################################################
### Model parametrization                                                      #
################################################################################

dataTrans <- function(data, x, y, z, tt, std.x, std.i, std.tt, inter, trace = TRUE){

  #########################################################
  ### Control and checking

  if(length(y) != 2)
    stop("\nThe outcome must contains 2 columns: 'time' and 'status'.")
  if(min(data[, y[1]]) < 0)
    stop("\nTimes must be non-negative.")
  if(sum(data[, y[2]] %in% c(0, 1)) != nrow(data))
    stop("\nStatus must be either 0 (censor) or 1 (event).")

  if(!is.null(tt)){
    TT <- as.data.frame(data[, tt])
    if(length(tt) > 1)
      stop("\nThe treatment must be a single variable.")
    if(length(tt) == 1){
      if(length(unique(as.vector(t(data[, tt])))) > 2)
        stop(paste0("\nThe treatment variable must consider only two groups."))
      if(length(unique(as.vector(t(data[, tt])))) == 1){
        warning(paste0("\nAll patients are in the same treatment group. The analysis is then switch to a prognostic setting."))
        tt <- NULL
        inter <- FALSE
      }
      if((sum(unique(as.vector(t(data[, tt]))) %in% c(-0.5, +0.5)) != 2) & std.tt == TRUE)
        data[, tt] <- as.numeric(factor(as.vector(t(data[, tt])))) - 1.5
    }
  }

  #########################################################

  itt <- tt; iz <- z; ix <- x; iy <- y
  iptt <- which(colnames(data) %in% tt)
  ipz <- which(colnames(data) %in% z)
  ipx <- which(colnames(data) %in% x)
  ipy <- which(colnames(data) %in% y)
  isSim <- (!is.null(attributes(data)$isSim))

  if(std.x == TRUE)
    data[, x] <- scale(data[, x], center = T, scale = T)

  if(inter == TRUE){
    XT <- as.matrix(data[, x]) * matrix(data[, tt], nrow = nrow(data), ncol = length(x))
    if(std.i == TRUE)
      XT <- scale(XT, center = T, scale = T)
    colnames(XT) <- paste0("bi", gsub(" ", "0", format(c(length(x), 1:length(x))))[-1])
    xt <- colnames(XT)
  }else{
    xt <- NULL
  }

  data <- cbind(data[, c(tt, z, x, y)])
  if(inter == TRUE){
    data <- cbind(data, XT)
    colnames(data)[1] <- "treat"
  }

  tnames <- c(rep("tt", length(tt)), rep("z", length(z)), rep("x", length(x)), rep("y", length(y)), rep("xt", length(xt)))
  if(!is.null(z))
  colnames(data)[tnames == "z"] <- paste0("cl", gsub(" ", "0", format(c(length(z), 1:length(z))))[-1])
  colnames(data)[tnames == "x"] <- paste0("bm", gsub(" ", "0", format(c(length(x), 1:length(x))))[-1])
  colnames(data)[tnames == "y"] <- c("time", "status")

  data <- na.omit(data)
  if(!is.null(attributes(data)$na.action) & trace == TRUE){
    nmiss <- length(attributes(data)$na.action)
    message(paste0(
      "\rData management: ", nmiss, " observation", ifelse(nmiss > 1, "s were", " was"), " excluded due to missing data."))
  }

  if(!(class(unlist(data)) %in% c("numeric", "integer")))
    stop("\nAll variables must be numerical.")

  attributes(data) <- append(
    x = attributes(data),
    values = list(
      inter = inter,
      inames = list(
        tt = itt,
        z = iz,
        x = ix,
        y = iy),
      ipos = list(
        tt = iptt,
        z = ipz,
        x = ipx,
        y = ipy
        ),
      tnames = tnames,
      pos = list(
        z = grep("cl", colnames(data)),
        x = grep("bm", colnames(data)),
        xt = grep("bi", colnames(data)),
        X = (1:ncol(data))[-which(colnames(data) %in% c("time", "status"))],
        y = which(colnames(data) %in% c("time", "status"))
        ),
      weights = c(rep(0, length(c(tt, z))), rep(1, length(x) * ((inter == TRUE) + 1))),
      std.x = std.x,
      std.tt = std.tt)
    )
    if(inter == TRUE){
      attributes(data)$pos$tt <- grep("treat", colnames(data))
      attributes(data)$std.i <- std.i
      attributes(data)$inames$xt <- paste0(ix, ":", itt)
    }

  return(data)
}
