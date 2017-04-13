################################################################################
################################################################################
################################################################################
### Prediction accuracy of the methods                                         #
################################################################################

predRes <- function(res, method, traindata, newdata, is.2cv, fold.2cv = 5, horizon, trace = TRUE, ncores = 1){

  ##########################################################
  ### Data checking and manipulation

  if(class(res) != "resBMsel")
    stop("\nThe 'res' object must be an object returned by the function BMsel().")

  if(missing(method)){
    method <- colnames(summary(res, show = FALSE, add.ridge = !is.na(attributes(res)$ridge)))
  }else{
    if(length(setdiff(method, names(res))) > 0)
      stop("\n Some methods in 'method' were not previously computed or do not exist.")
    method <- unique(c(method, if(attributes(res)$isSim == TRUE) "oracle"))
  }

  if(missing(traindata))
    stop("\nThe training data set used in the BMsel() function must be specified.")
  traindata <- as.data.frame(traindata)
  rownames(traindata) <- 1:nrow(traindata)

  if(length(setdiff(attributes(res)$inames[which(attributes(res)$tnames != 'xt')], colnames(traindata))) > 0)
    stop("\nSome covariates of the training set are missing. Please specify the same data set used for BMsel().")

  if(missing(is.2cv))
    is.2cv <- FALSE

  if(!(is.2cv) %in% c(TRUE, FALSE))
    stop("\nThe 'is.2cv' parameter must be either TRUE or FALSE.")

  if(fold.2cv < 2 || fold.2cv > nrow(traindata))
    stop("\nThe 'fold.2cv' parameter must be between 2 and the sample size 'n'.")

  if(missing(horizon)){
    stop("\nThe 'horizon' parameter must be specified for Cox models.")
  }else{
    if(min(horizon) < 0 || max(horizon) > min(c(max(traindata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')[1]]]),
                                        if(!missing(newdata)) max(newdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')[1]]]))))
      stop("\n The 'horizon' parameter is out of the range of the observed survival time.")
  }
  horizon <- sort(horizon)

  if(!missing(newdata)){
    if(length(setdiff(attributes(res)$inames[which(attributes(res)$tnames != 'xt')], colnames(newdata))) > 0)
      stop("Some covariates of the new data set are missing. Please specify the same covariates as the training set.")
    newdata <- as.data.frame(newdata)
    rownames(newdata) <- 1:nrow(newdata)
  }

  ncores <- round(ncores, 0)
  if(ncores < 1 || ncores > detectCores())
    stop(paste0("\n'ncores' must be between 1 and ", detectCores(), "."))
  if(ncores > fold.2cv) ncores <- fold.2cv

  ##########################################################
  ### Data management

  tt <- attributes(res)$inames[which(attributes(res)$tnames == 'tt')]
  x <- attributes(res)$inames[which(attributes(res)$tnames == 'x')]
  z <- attributes(res)$inames[which(attributes(res)$tnames == 'z')]
  y <- attributes(res)$inames[which(attributes(res)$tnames == 'y')]
  isRidge <- (unique(!is.na(attributes(res)$ridge)))
  isNew <- (!missing(newdata))

  Res <- data.frame(summary(res, show = FALSE, add.ridge = isRidge))
  Res <- Res[, gsub("-", ".", method), drop = FALSE]

  if(attributes(res)$inter == TRUE){
    Res.i <- data.frame(summary(res, show = FALSE, keep = "xt", add.ridge = isRidge))
    Res.i <- Res.i[, gsub("-", ".", method), drop = FALSE]
  }

  nmeth <- ncol(Res)

  ##########################################################
  ### Training set
  if(trace == TRUE)
    message(paste0(
      "\rComputing prediction criteria for: training set"))

  attr.train <- attributes(traindata)[-which(names(attributes(traindata)) %in% "names")]

  tdata <- dataTrans(data = traindata, x = x, y = y, z = z, tt = tt,
                     std.x = attributes(res)$std.x, std.i = attributes(res)$std.i, std.tt = attributes(res)$std.tt,
                     inter = attributes(res)$inter, trace = FALSE)
  colnames(tdata) <- attributes(res)$inames
  attributes(tdata) <- append(attributes(tdata), attr.train)

  surv.train <- Surv(tdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')][1]],
                     tdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')][2]])

  lp.train <- matrix(0, nrow = nrow(tdata), ncol = ncol(Res))
  if(nrow(Res) > 0)
    lp.train <- as.matrix(tdata[, rownames(Res)]) %*% as.matrix(Res)

  if(attributes(res)$inter == TRUE){
    lpint.train <- matrix(0, nrow = nrow(tdata), ncol = ncol(Res))
    if(nrow(Res.i) > 0 & sum(Res.i) != 0)
      lpint.train <- as.matrix(tdata[, gsub(".INT", "", gsub("bi", "bm", rownames(Res.i)))]) %*% as.matrix(Res.i)
  }

  predRes.train <- compute.predRes(res = res, nmeth = nmeth, hrz = horizon, traindata = traindata, newdata = traindata,
                                   surv.train = surv.train, surv.new = surv.train, lp.train = lp.train, lp.new = lp.train,
                                   lpint.train = lpint.train, lpint.new = lpint.train, tt = tt)

  ##########################################################
  ### External validation
  if(!missing(newdata)){
    if(trace == TRUE)
      message(paste0(
        "\rComputing prediction criteria for: validation set"))

    newdata <- dataTrans(data = newdata, x = x, y = y, z = z, tt = tt,
                         std.x = attributes(res)$std.x, std.i = attributes(res)$std.i, std.tt = attributes(res)$std.tt,
                         inter = attributes(res)$inter, trace = TRUE)
    colnames(newdata) <- attributes(res)$inames

    surv.new <- Surv(newdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')][1]],
                     newdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')][2]])

    lp.new <- matrix(0, nrow = nrow(newdata), ncol = ncol(Res))
    if(nrow(Res) > 0)
      lp.new <- as.matrix(newdata[, rownames(Res)]) %*% as.matrix(Res)

    if(attributes(res)$inter == TRUE){
      lpint.new <- matrix(0, nrow = nrow(newdata), ncol = ncol(Res))
      if(nrow(Res.i) > 0 & sum(Res.i) != 0)
        lpint.new <- as.matrix(newdata[, gsub(".INT", "", gsub("bi", "bm", rownames(Res.i)))]) %*% as.matrix(Res.i)
    }

    predRes.new <- compute.predRes(res = res, nmeth = nmeth, hrz = horizon, traindata = traindata, newdata = newdata,
                                   surv.train = surv.train, surv.new = surv.new, lp.train = lp.train, lp.new = lp.new,
                                   lpint.train = lpint.train, lpint.new = lpint.new, tt = tt)
  }

  ##########################################################
  ### Internal evaluation through double cross-validation (2cv)
  if(is.2cv == TRUE){
    form <- attributes(res)$formula
    foldid2 <- sample(x = 1:fold.2cv, size = nrow(traindata), replace = T)

    if(trace == TRUE)
      message(
        "\rComputing prediction criteria for: internal 2CV")

    cl <- makeCluster(ncores)

    res.2cv <- clusterApplyLB(
      cl = cl,
      x = 1:fold.2cv,
      fun = function(X){

        traindataT <- traindata[which(foldid2 != X), ]
        ltraindataT <- list(traindataT)
        traindataV <- traindata[which(foldid2 == X), ]

        form[2] <- ltraindataT
        form[which(names(form) == "method")] <- list(setdiff(method, "oracle"))
        w <- options()$warn
        options(warn = -1)
        pos.trace <- which(names(form) == "trace")
        if(length(pos.trace) == 0){
          form[length(names(form)) + 1] <- FALSE
          names(form)[length(names(form))] <- "trace"
        }else{
          form[pos.trace] <- FALSE
        }
        res <- eval(form)
        options(warn = w)

        Res <- data.frame(summary(res, show = FALSE, add.ridge = isRidge))
        Res <- Res[, gsub("-", ".", method), drop = FALSE]

        if(attributes(res)$inter == TRUE){
          Res.i <- data.frame(summary(res, show = FALSE, keep = "xt", add.ridge = isRidge))
          Res.i <- Res.i[, gsub("-", ".", method), drop = FALSE]
        }

        nmeth <- ncol(Res)

        traindataV <- dataTrans(data = traindataV, x = x, y = y, z = z, tt = tt,
                                std.x = attributes(res)$std.x, std.i = attributes(res)$std.i, std.tt = attributes(res)$std.tt,
                                inter = attributes(res)$inter, trace = FALSE)
        colnames(traindataV) <- attributes(res)$inames

        lp.trainV <- matrix(0, nrow = nrow(traindataV), ncol = ncol(Res))
        if(nrow(Res) > 0)
          lp.trainV <- as.matrix(traindataV[, rownames(Res)]) %*% as.matrix(Res)

        lpint.trainV <- NA
        if(attributes(res)$inter == TRUE){
          lpint.trainV <- matrix(0, nrow = nrow(traindataV), ncol = ncol(Res))
          if(nrow(Res.i) > 0 & sum(Res.i) != 0)
            lpint.trainV <- as.matrix(traindataV[, gsub(".INT", "", gsub("bi", "bm", rownames(Res.i)))]) %*% as.matrix(Res.i)
        }

        rownames(lp.trainV) <- rownames(traindataV)
        colnames(lp.trainV) <- method

        if(attributes(res)$inter == TRUE){
          rownames(lpint.trainV) <- rownames(traindataV)
          colnames(lpint.trainV) <- method
        }

        return(list(lp.trainV, lpint.trainV))
      }
    )

    stopCluster(cl)

    lp.2cv <- data.frame()
    for(i in 1:fold.2cv)
      lp.2cv <- rbind(lp.2cv, data.frame(res.2cv[[i]][1]))
    lp.2cv <- lp.2cv[order(as.numeric(rownames(lp.2cv))), , drop = FALSE]

    if(attributes(res)$inter == TRUE){
      lpint.2cv <- data.frame()
      for(i in 1:fold.2cv)
        lpint.2cv <- rbind(lpint.2cv, data.frame(res.2cv[[i]][2]))
      lpint.2cv <- lpint.2cv[order(as.numeric(rownames(lpint.2cv))), , drop = FALSE]
    }

    predRes.2cv <- compute.predRes(res = res, nmeth = nmeth, hrz = horizon, traindata = traindata, newdata = traindata,
                                   surv.train = surv.train, surv.new = surv.train, lp.train = lp.train, lp.new = lp.2cv,
                                   lpint.train = lpint.train, lpint.new = lpint.2cv, tt = tt)

  }

  ##########################################################
  predRes <- list()
  for(i in 1:length(horizon)){
    predres <- list('Training set' = round(predRes.train[[i]], 4))
    if(is.2cv == TRUE) predres <- merge.list(predres, list('Internal validation (2CV)' = round(predRes.2cv[[i]], 4)))
    if(!missing(newdata)) predres <- merge.list(predres, list('External validation' = round(predRes.new[[i]], 4)))
    predRes[[i]] <- predres
  }
  names(predRes) <- paste0("horizon = ", horizon)
  class(predRes) <- "predRes"
  return(predRes)
}

compute.predRes <- function(res, nmeth, hrz, traindata, newdata, surv.train, surv.new,
                 lp.train, lp.new, lpint.train, lpint.new, tt){
  nres <- colnames(lp.train)
  pRes <- lapply(
    X = 1:length(hrz),
    FUN = function(X){
      hrz <- hrz[X]
      pres <- matrix(unlist(lapply(
        X = 1:nmeth,
        FUN = function(X){
          if(sum(lp.new[, X] != 0) > 0){
            uno <- UnoC(
              Surv.rsp = surv.train,
              Surv.rsp.new = surv.new,
              lpnew = lp.new[, X] + rnorm(length(lp.new[, X]), 0, 1e-7),
              time = hrz)
          }else{
            uno <- 0.5
          }
          brier <- predErr(
            Surv.rsp = surv.train,
            Surv.rsp.new = surv.new,
            lp = lp.train[, X],
            lpnew = lp.new[, X],
            times = seq(0, hrz, hrz/20),
            type = "brier",
            int.type = "weighted")$i
          if(attributes(res)$inter == TRUE){
            if(sum(lpint.new[, X] != 0) > 0){
              pC <- which(traindata[, tt] == min(unique(traindata[, tt]))); pC.New <- which(newdata[, tt] == min(unique(newdata[, tt])))
              pT <- which(traindata[, tt] == max(unique(traindata[, tt]))); pT.New <- which(newdata[, tt] == max(unique(newdata[, tt])))
              duno <- abs(
                UnoC(
                  Surv.rsp = surv.train[pT],
                  Surv.rsp.new = surv.new[pT.New],
                  lpnew = lpint.new[pT.New, X] + rnorm(length(lpint.new[pT.New, X]), 0, 1e-7),
                  time = hrz) -
                  UnoC(
                    Surv.rsp = surv.train[pC],
                    Surv.rsp.new = surv.new[pC.New],
                    lpnew = lpint.new[pC.New, X] + rnorm(length(lpint.new[pC.New, X]), 0, 1e-7),
                    time = hrz))
            }else{
              duno <- 0
            }
          }
          return(c(uno, brier, if(attributes(res)$inter == TRUE) duno))
        })), nrow = ifelse(attributes(res)$inter == TRUE, 3, 2), ncol = nmeth, byrow = FALSE)
      colnames(pres) <- gsub("[.]", "-", nres)
      rownames(pres) <- c("C-index", "Prediction Error", if(attributes(res)$inter == TRUE) "Delta C-index")
      return(pres)
    })
  names(pRes) <- hrz
  return(pRes)
}
###########################################################

print.predRes <- function(x, ...) {

  n <- names(x)
  attributes(x) <- NULL
  names(x) <- n
  print(x, ...)

}

###########################################################

plot.predRes <- function(x, method, crit = c("C", "PE", "dC"), xlim, ylim, xlab, ylab, col, ...){

  if(missing(method)){
    method <- colnames(x[[1]][[1]])
  }else{
    if(length(setdiff(method, colnames(x[[1]][[1]]))) > 0)
      stop("Some methods in 'method' do not exist.")
  }

  if(missing(crit)){
    stop("'crit' must be specified.")
  }else{
    crit <- match.arg(crit, several.ok = FALSE)
    if(crit == "dC" & nrow(x[[1]][[1]]) == 2)
      stop("'crit = dC' can not be specified for a prognostic setting.")
  }
  crit = switch(crit, C = "C-index", PE = "Prediction Error", dC = "Delta C-index")

  if(missing(col)){
    col <- 1:length(method)
  }else{
    if(!length(col) %in% c(1, length(method)))
      stop("'col' must be of length 1 or the number of methods 'method'.")
  }

  hrz <- as.numeric(gsub("horizon = ", "", names(x)))

  xx <- data.frame()
  for(i in 1:length(x)){
    for(j in 1:length(x[[i]])){
      pos <- (j-1) * length(method) + 1:length(method)
      xx[pos, 1] <- switch(names(x[[i]][j]), 'Training set' = 3, 'Internal validation (2CV)' = 2, 'External validation' = 1)
      xx[pos, 2] <- method
      xx[pos, 2 + i] <- x[[i]][[j]][crit, method]
    }
  }
  xx[, 2] <- col[match(xx[, 2], method)]

  if(missing(xlim))
    xlim <- range(hrz)

  if(missing(ylim)){
    if(crit == "C-index"){
      ylim <- c(0.5, 1)
    }else{
      ylim <- c(0, 0.5)
    }
  }

  split.screen(rbind(c(0, 0.33, 0.85, 1),
                     c(0.33, 1, 0.85, 1),
                     c(0, 1, 0, 0.85)
  ))

  screen(1)
  par(mar = rep(0, 4))
  frame()
  legend(0, .5, xjust = 0, yjust = .5, legend = names(x[[1]]), pch = NA, lty = unique(xx[, 1]), ncol = 1, box.col = "white", pt.cex = 1.2, cex = 1.2)

  screen(2)
  par(mar = rep(0, 4))
  frame()
  legend(0, .5, xjust = 0, yjust = .5, legend = method, pch = rep(22, 5), lty = rep(0, length(method)), ncol = 5, col = unique(xx[, 2]), pt.bg = unique(xx[, 2]), box.col = "white", pt.cex = 1.2, cex = 1.2)

  screen(3)
  par(mar = c(5, 5, 0, 3))
  plot(NULL, xlim = xlim, ylim = ylim, xlab = "Time", ylab = crit, ...)
  zz <- lapply(
    X = 1:nrow(xx),
    FUN = function(X){
      lines(x = hrz, y = xx[X, -(1:2)], lty = xx[X, 1], col = xx[X, 2], ...)
    })
  close.screen(all.screens = TRUE)

} # end of plot
###########################################################
