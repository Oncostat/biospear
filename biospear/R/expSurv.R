################################################################################
################################################################################
################################################################################
### Global function for estimated expected survival probabilities              #
################################################################################

expSurv <- function(
  ####################################################################
  ######################## *** PARAMETERS *** ########################
  res,              # Object resulting from the BMsel() function
  traindata,        # Data set used to compute the 'res' object from BMsel()
  method,           # Methods to compute
  ci.level = .95,   # Level for a two-sided confidence interval
  boot = FALSE,     # Computation of confidence intervals via bootstrap? (logical: TRUE/FALSE)
  nboot,            # Number of bootstrap samples
  smooth,           # Computation of smoothed B-splines? (logical: TRUE/FALSE)
  pct.group,        # Number or percentile of the prognostic category (interaction setting only)
  horizon,          # Time-horizon to estimate the expected survival probabilities
  trace = TRUE,     # Print messages ?
  ncores = 1        # Number of PC cores used (useful for the bootstrap CI only)
  #####################################################################
  ) {

  #####################################################################
  ### CONTROL CHECKING

  if(length(setdiff(attributes(res)$inames[which(attributes(res)$tnames != 'xt')], colnames(traindata))) > 0)
    stop("\nSome covariates of the 'traindata' set are missing. Please specify the training set used in BMsel().")

  if(missing(method)){
    method <- colnames(summary(res, show = FALSE, add.ridge = !is.na(attributes(res)$ridge)))
  }else{
    if(length(setdiff(method, names(res))) > 0)
      stop("\n Some methods in 'method' were not previously computed or do not exist.")
    method <- unique(c(method, if(attributes(res)$isSim == TRUE) "oracle"))
  }

  if(!is.logical(boot) || length(boot) > 1)
    stop("\n'boot' must be TRUE or FALSE.")

  if(missing(nboot)){
    if(boot == TRUE)
      stop("\n'nboot' must be specified (as 'boot' = TRUE).")
  }else{
    if(boot == FALSE)
      warning("\n'nboot' is not considered (as 'boot' = FALSE).")
    if(nboot < 0 || !is.numeric(nboot) || length(boot) > 1)
      stop("\n'nboot' must be a positive integer.")
    nboot <- ceiling(nboot)
  }

  if(!is.logical(smooth) || length(boot) > 1)
    stop("\n'smooth' must be TRUE or FALSE.")

  if(attributes(res)$inter == FALSE){
    if(!missing(pct.group))
      warning("\nThe 'pct.group' parameter is not considered for the prognostic setting.")
    pct.group <- 1
  }else{
    if(missing(pct.group))
      stop("\n'pct.group' must be specified for the interaction setting.")
    if(!is.numeric(pct.group))
      stop("\n'pct.group' must be numerical.")
    if(length(pct.group) == 1){
      if(pct.group < 0 || pct.group > 6)
        stop("\nThe number of prognostic groups must be between 1 and 6.")
      pct.group <- ceiling(pct.group)
      pct.group <- switch(
        pct.group,
        "1" = c(0, 1),
        "2" = c(0, 0.5, 1),
        "3" = c(0, 0.27, 0.73, 1),
        "4" = c(0, 0.164, 0.5, 0.886, 1),
        "5" = c(0, 0.11, 0.346, 0.653, 0.89, 1),
        "6" = c(0, 0.074, 0.255, 0.5, 0.755, 0.926, 1)
      )
    }else{
      if(sum(pct.group) != 1)
        stop("The sum of the prognostic group percentiles in 'pct.group' is not equal to 1.")
      pct.group <- c(0, cumsum(pct.group))
    }
  }

  if(missing(horizon)){
    stop("\'horizon' must be specified for Cox models.")
  }else{
    if(length(horizon) > 1)
      stop("\n'horizon' must be an unique value.")
    if(horizon < 0 || horizon > max(traindata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')[1]]]))
      stop("\n'horizon' is out of the range of the observed survival time.")
  }

  if(!is.logical(trace))
    stop("\n'trace' must be either TRUE or FALSE.")

  ncores <- round(ncores, 0)
  if(ncores < 1 || ncores > detectCores())
    stop(paste0("\n'ncores' must be between 1 and ", detectCores(), "."))

  #####################################################################
  ### COMPUTATION OF THE PREDICTION SCORES

  isRidge <- unique(!is.na(attributes(res)$ridge))
  m.sel <- data.frame(summary(res, show = FALSE, keep = c("t", "z", "x"), add.ridge = isRidge))
  colnames(m.sel) <- gsub("[.]", "-", colnames(m.sel))
  m.sel <- m.sel[, method, drop = FALSE]

  m.score <- matrix(0, nrow = nrow(traindata), ncol = ncol(m.sel))
  if(nrow(m.sel) > 0 & sum(m.sel) != 0){
    m.score <- as.matrix(traindata[, rownames(m.sel), drop = TRUE]) %*% as.matrix(m.sel)
  }
  m.score[is.na(m.score)] <- 0

  if(attributes(res)$inter == TRUE){

    i.sel <- data.frame(summary(res, show = FALSE, keep = "xt", add.ridge = isRidge))
    colnames(i.sel) <- gsub("[.]", "-", colnames(i.sel))
    i.sel <- i.sel[, method, drop = FALSE]
    i.score <- matrix(0, nrow = nrow(traindata), ncol = ncol(i.sel))

    if(nrow(i.sel) > 0 & sum(i.sel) != 0){

      i.score <- as.matrix(traindata[, gsub(paste0(":", attributes(res)$inames[1]), "", rownames(i.sel))]) %*% as.matrix(i.sel)
      i.scale <- apply(i.score, 2, FUN = function(X) quantile(X, probs = c(0.025, 0.975)))
      i.score <- matrix(unlist(lapply(
        X = 1:ncol(i.scale),
        FUN = function(X){
          (i.score[, X] - i.scale[1, X]) / diff(i.scale[, X])
      })), nrow = nrow(traindata), ncol = ncol(i.sel))
      colnames(i.score) <- colnames(i.sel)
      i.score[is.na(i.score)] <- 0

      m.scale <- NULL
      m.breaks <- matrix(unlist(lapply(
        X = 1:ncol(m.score),
        FUN = function(X){
          if(length(pct.group) - 1 > length(unique(m.score[, X]))){
            warning(paste0(
              "\nThe number of prognostic group (",
              length(pct.group) - 1,
              ") is larger than the number of values for '",
              colnames(m.score)[X],
              "' (", length(unique(m.score[, X])), ")"))
            c(sort(unique(m.score[, X])), rep(NA, length(pct.group) - length(unique(m.score[, X]))))
          }else{
            quantile(m.score[, X], probs = pct.group)
          }
        })), nrow = length(pct.group), byrow = FALSE)
      colnames(m.breaks) <- method
      catm.score <- matrix(unlist(lapply(
        X = 1:ncol(m.score),
        FUN = function(X){
          if(sum(m.score[, X] != 0) > 0){
            if(sum(is.na(m.breaks[, X])) > 0){
              match(m.score[, X], na.exclude(m.breaks[, X]))
            }else{
              cut(m.score[, X], breaks = m.breaks[, X], labels = FALSE, include.lowest = TRUE)
            }
          }else{
            rep(1, length(m.score[, X]))
          }
        }
      )), ncol = ncol(m.sel), byrow = FALSE, dimnames = list(1:nrow(traindata), colnames(m.sel)))

    }else{
      warning("\nAs no interaction was selected, we fixed 'pct.group' = 1 and 'inter' = FALSE.")
      pct.group <- 1
      attributes(res)$inter <- FALSE
    }
  }
  if(attributes(res)$inter == FALSE){
    i.scale <- i.score <- i.sel <- m.breaks <- catm.score <- NULL
    m.scale <- apply(m.score, 2, FUN = function(X) quantile(X, probs = c(0.025, 0.975)))
    m.score <- matrix(unlist(lapply(
      X = 1:ncol(m.scale),
      FUN = function(X)
        (m.score[, X] - m.scale[1, X]) / diff(m.scale[, X])
      )), nrow = nrow(traindata), ncol = ncol(m.sel))
    colnames(m.score) <- colnames(m.sel)
    m.score[is.na(m.score)] <- 0
  }
  #####################################################################
  ### COMPUTATION OF THE EXPECTED SURVIVAL
  ### + ANALYTICAL CI, if 'boot' = FALSE

  all.sel <- data.frame(summary(res, show = FALSE, add.ridge = unique(!is.na(attributes(res)$ridge))))
  colnames(all.sel) <- gsub("[.]", "-", colnames(all.sel))
  all.sel <- all.sel[, method, drop = FALSE]

  if(trace == TRUE){
    message("\rComputation of the expected survival")
    if(boot == FALSE)
      message("\rComputation of analytical confidence intervals")
    flush.console()
  }

  fct <- as.formula(paste0("Surv(", paste(
    attributes(res)$inames[attributes(res)$tnames == "y"], collapse = ", "),
    ") ~ 1"))
  sbhz <- survfit(fct, data = traindata, se.fit = TRUE, conf.int = ci.level)
  bhz <- c(
    surv = sbhz$surv[max(which(sbhz$time < horizon))],
    lower = sbhz$lower[max(which(sbhz$time < horizon))],
    upper = sbhz$upper[max(which(sbhz$time < horizon))])

  surv <- int.expSurv(res = res, all.sel = all.sel, datafit = traindata, dataest = traindata,
                           boot = boot, ci.level = ci.level, hrz = horizon, out.fit = TRUE)
  m <- matrix(0, nrow = nrow(traindata), ncol = ncol(all.sel))
  colnames(m) <- colnames(all.sel)
  train.surv <- list(surv = m, lower = m, upper = m)
  train.fit <- list()
  for(i in 1:length(surv)){
    train.surv$surv[, i] <- surv[[i]][[1]]$surv
    if(boot == FALSE) train.surv$lower[, i] <- surv[[i]][[1]]$lower
    if(boot == FALSE) train.surv$upper[, i] <- surv[[i]][[1]]$upper
    train.fit[[i]] <- surv[[i]][[2]]
  }

  #####################################################################
  ### COMPUTATION OF BOOTSTRAPED CI (if 'boot' = TRUE)

  mat.boot <- NULL
  if(boot == TRUE){
    form <- attributes(res)$formula
    if(trace == TRUE)
      message("\rComputation of bootstrap confidence intervals")

    cl <- makeCluster(ncores)

    surv.boot <- clusterApplyLB(
      cl = cl,
      x = 1:nboot,
      fun = function(X){

        if(trace == TRUE)
          message(paste0("\rBootstrap sample ", X, " / ", nboot))

        bootdata <- traindata[sample(1:nrow(traindata), size = nrow(traindata), replace = TRUE), ]
        form[2] <- list(bootdata)
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

        all.sel <- data.frame(summary(res, show = FALSE, add.ridge = unique(!is.na(attributes(res)$ridge))))
        colnames(all.sel) <- gsub("[.]", "-", colnames(all.sel))
        all.sel <- all.sel[, method, drop = FALSE]

        surv <- int.expSurv(res = res, all.sel = all.sel, datafit = bootdata, dataest = traindata,
                            boot = boot, ci.level = ci.level, hrz = horizon, out.fit = TRUE)

        for(i in 1:ncol(all.sel)){
          boot.bhz <- basehaz(surv[[i]][[2]], centered = T)
          surv[[i]][[2]] <- c(surv.bhz = boot.bhz[max(which(boot.bhz[, "time"] < horizon)), "hazard"], coef(surv[[i]][[2]]))
        }

        return(surv)
      })

    stopCluster(cl)

    mat.boot <- list()
    for(i in 1:ncol(all.sel)){
      trainb.mat <- matrix(unlist(lapply(
        X = 1:nboot,
        FUN = function(X){
            return(surv.boot[[X]][[i]][[1]]$surv)
          })), nrow = nrow(traindata), ncol = nboot, byrow = FALSE)
      train.surv$lower[, i] <- apply(X = trainb.mat, MARGIN = 1, FUN = function(X) quantile(X, (1 - ci.level) / 2))
      train.surv$upper[, i] <- apply(X = trainb.mat, MARGIN = 1, FUN = function(X) quantile(X, 1 - (1 - ci.level) / 2))

      names <- unique(names(unlist(sapply(surv.boot, "[[", i)[2, ])))
      m.boot <- matrix(0, ncol = length(names), nrow = nboot, dimnames = list(1:nboot, names))
      for(X in 1:nboot){
        r <- sapply(surv.boot, "[[", i)[2, ][[X]]
        m.boot[X, names(r)] <- r
      }
      mat.boot[[i]] <- m.boot
    }

  }

  #####################################################################
  ## COMPUTATION OF SMOOTHED B-SPLINES

  strain.surv <- mat.smooth <- NULL

  if(smooth == TRUE){

    if(trace == TRUE)
      message("\rComputation of smoothed B-splines")

    strain.surv <- train.surv
    nmeth <- ncol(m.score)

    if(attributes(res)$inter == FALSE){

      mat.smooth <- matrix(0, nrow = 3 * nmeth, ncol = 13,
                           dimnames = list(paste0(rep(colnames(m.score), each = 3), c(".surv", ".lower", ".upper")),
                                           c(paste0("knots", 1:6), paste0("coef", 1:7))))
      for(i in 1:nmeth){

        if(length(unique(m.score[, i])) > 2){
          sres <- int.expSSurv(x = m.score[, i], y = train.surv$lower[, i], z = m.score[, i], c = "decrease")
          strain.surv$lower[, i] <- sres[[1]]; mat.smooth[3 * (i-1) + 2, names(sres[[2]])] <- sres[[2]]
          sres <- int.expSSurv(x = m.score[, i], y = train.surv$upper[, i], z = m.score[, i], c = "decrease")
          strain.surv$upper[, i] <- sres[[1]]; mat.smooth[3 * (i-1) + 3, names(sres[[2]])] <- sres[[2]]
        }

      }

    }else{

      ncat <- (length(pct.group) - 1)
      mat.smooth <- matrix(0, nrow = 6 * nmeth * ncat, ncol = 13,
                           dimnames = list(paste0(rep(colnames(m.score), each = 6 * ncat), ".",
                                                  rep(1:ncat, each = 6),
                                                  rep(c("t", "c"), each = 3), ".",
                                                  c("surv", "lower", "upper")),
                                           c(paste0("knots", 1:6), paste0("coef", 1:7))))

      for(i in 1:nmeth){
        for(j in 1:ncat){

          wht <- which(catm.score[, i] == j & traindata[, attributes(res)$inames[1]] == +0.5)
          whc <- which(catm.score[, i] == j & traindata[, attributes(res)$inames[1]] == -0.5)
          pos <- paste0(colnames(m.score)[i], ".", j)

          if(sum(i.score[wht, i] != 0) > 0 & sum(i.score[whc, i] != 0) > 0){

            sres <- int.expSSurv(x = i.score[wht, i], y = train.surv$surv[wht, i], z = i.score[wht, i], c = "decrease")
            strain.surv$surv[wht, i] <- sres[[1]]; mat.smooth[paste0(pos, "t.surv"), names(sres[[2]])] <- sres[[2]]
            sres <- int.expSSurv(x = i.score[wht, i], y = train.surv$lower[wht, i], z = i.score[wht, i], c = "decrease")
            strain.surv$lower[wht, i] <- sres[[1]]; mat.smooth[paste0(pos, "t.lower"), names(sres[[2]])] <- sres[[2]]
            sres <- int.expSSurv(x = i.score[wht, i], y = train.surv$upper[wht, i], z = i.score[wht, i], c = "decrease")
            strain.surv$upper[wht, i] <- sres[[1]]; mat.smooth[paste0(pos, "t.upper"), names(sres[[2]])] <- sres[[2]]

            sres <- int.expSSurv(x = i.score[whc, i], y = train.surv$surv[whc, i], z = i.score[whc, i], c = "increase")
            strain.surv$surv[whc, i] <- sres[[1]]; mat.smooth[paste0(pos, "c.surv"), names(sres[[2]])] <- sres[[2]]
            sres <- int.expSSurv(x = i.score[whc, i], y = train.surv$lower[whc, i], z = i.score[whc, i], c = "increase")
            strain.surv$lower[whc, i] <- sres[[1]]; mat.smooth[paste0(pos, "c.lower"), names(sres[[2]])] <- sres[[2]]
            sres <- int.expSSurv(x = i.score[whc, i], y = train.surv$upper[whc, i], z = i.score[whc, i], c = "increase")
            strain.surv$upper[whc, i] <- sres[[1]]; mat.smooth[paste0(pos, "c.upper"), names(sres[[2]])] <- sres[[2]]
          }
        }
      }
    }
  }
  if(smooth == TRUE) surv <- strain.surv else surv <- train.surv

  tt <- NULL
  if(attributes(res)$tnames[1] == "tt")
    tt = traindata[, attributes(res)$inames[1]]

  attributes(surv) <- append(
    x = attributes(train.surv),
    values = list(
      all.sel = all.sel,
      m.sel = m.sel,
      m.breaks = m.breaks,
      m.score = m.score,
      catm.score = catm.score,
      m.scale = m.scale,
      i.sel = i.sel,
      i.score = i.score,
      i.scale = i.scale,
      tt = tt,
      inter = attributes(res)$inter,
      inames = attributes(res)$inames,
      tnames = attributes(res)$tnames,
      std.x = attributes(res)$std.x,
      std.i = attributes(res)$std.i,
      ci.level = ci.level,
      train.fit = train.fit,
      boot = boot,
      mat.boot = mat.boot,
      smooth = smooth,
      mat.smooth = mat.smooth,
      hrz = horizon,
      bhz = bhz
    ))
  if(smooth == TRUE)
    attributes(surv)$surv <- train.surv$surv

  class(surv) <- 'resexpSurv'
  return(surv)
}
###########################################################
###########################################################

print.resexpSurv <- function(x, ...){

  n <- names(x)
  attributes(x) <- NULL
  names(x) <- n
  print(x, ...)

}
###########################################################

predict.resexpSurv <- function(object, newdata, ...){

  a <- attributes(object)
  m <- matrix(0, nrow = nrow(newdata), ncol = ncol(a$all.sel))
  colnames(m) <- colnames(a$all.sel)
  new.surv <- list(surv = m, lower = m, upper = m)

  m.score <- matrix(0, nrow = nrow(newdata), ncol = ncol(a$m.sel))
  if(nrow(a$m.sel) > 0){
    m.score <- as.matrix(newdata[, rownames(a$m.sel), drop = TRUE]) %*% as.matrix(a$m.sel)
  }
  m.score[is.na(m.score)] <- 0
  if(a$inter == TRUE){
    catm.score <- matrix(unlist(lapply(
      X = 1:ncol(m.score),
      FUN = function(X)
        cut(m.score[, X], breaks = a$m.breaks[, X], labels = FALSE, include.lowest = TRUE)
    )), ncol = ncol(a$m.sel), byrow = FALSE, dimnames = list(1:nrow(newdata), colnames(a$m.sel)))
    i.score <- matrix(0, nrow = nrow(newdata), ncol = ncol(a$i.sel))
    if(nrow(a$i.sel) > 0){
      i.score <- as.matrix(newdata[, gsub(paste0(":", attributes(object)$inames[1]), "", rownames(a$i.sel))]) %*% as.matrix(a$i.sel)
    }
    i.score <- matrix(unlist(lapply(
      X = 1:ncol(a$i.scale),
      FUN = function(X)
        (i.score[, X] - a$i.scale[1, X]) / diff(a$i.scale[, X])
    )), nrow = nrow(newdata), ncol = ncol(a$i.scale))
    colnames(i.score) <- colnames(a$i.scale)
  }else{
    i.score <- catm.score <- NULL
    m.score <- matrix(unlist(lapply(
      X = 1:ncol(a$m.scale),
      FUN = function(X)
        (m.score[, X] - a$m.scale[1, X]) / diff(a$m.scale[, X])
    )), nrow = nrow(newdata), ncol = ncol(a$m.scale))
    colnames(m.score) <- colnames(a$m.scale)
    m.score[is.na(m.score)] <- 0
  }

  for(i in 1:ncol(a$all.sel)){
    if(a$smooth == TRUE){
      if(attributes(object)$inter == TRUE){
        for(j in 1:max(catm.score[, i], na.rm = T)){
          wht <- which(catm.score[, i] == j & newdata[, a$inames[1]] == +0.5)
          whc <- which(catm.score[, i] == j & newdata[, a$inames[1]] == -0.5)
          pos <- paste0(colnames(m.score)[i], ".", j)

          if(length(unique(m.score[, i])) > 2){

            new.surv$surv[wht, i] <- int.predSmooth(
              i = i, label = paste0(pos, "t.surv"), score = i.score, wrow = wht, mat = a$mat.smooth)
            new.surv$surv[whc, i] <- int.predSmooth(
              i = i, label = paste0(pos, "c.surv"), score = i.score, wrow = whc, mat = a$mat.smooth)

            new.surv$lower[wht, i] <- int.predSmooth(
              i = i, label = paste0(pos, "t.lower"), score = i.score, wrow = wht, mat = a$mat.smooth)
            new.surv$lower[whc, i] <- int.predSmooth(
              i = i, label = paste0(pos, "c.lower"), score = i.score, wrow = whc, mat = a$mat.smooth)

            new.surv$upper[wht, i] <- int.predSmooth(
              i = i, label = paste0(pos, "t.upper"), score = i.score, wrow = wht, mat = a$mat.smooth)
            new.surv$upper[whc, i] <- int.predSmooth(
              i = i, label = paste0(pos, "c.upper"), score = i.score, wrow = whc, mat = a$mat.smooth)
          }
        }
      }else{

        if(sum(m.score[, i] != 0) > 0){
          if(!(length(unique(m.score[, i])) > 2)){
            s <- survfit(a$train.fit[[i]], newdata = newdata, se.fit = !(length(unique(m.score[, i])) > 2))
            new.surv$surv[, i] <- s$surv[max(which(s$time < a$hrz)), ]
            new.surv$lower[, i] <- s$lower[max(which(s$time < a$hrz)), ]
            new.surv$upper[, i] <- s$upper[max(which(s$time < a$hrz)), ]
          }else{
            pos <- colnames(m.score)[i]
            new.surv$lower[, i] <- int.predSmooth(
              i = i, label = paste0(pos, ".lower"), score = m.score, mat = a$mat.smooth)
            new.surv$upper[, i] <- int.predSmooth(
              i = i, label = paste0(pos, ".upper"), score = m.score, mat = a$mat.smooth)
          }
        }else{
          new.surv$surv[, i] <- rep(a$bhz["surv"], nrow(newdata))
          new.surv$lower[, i] <- rep(a$bhz["lower"], nrow(newdata))
          new.surv$upper[, i] <- rep(a$bhz["upper"], nrow(newdata))
        }
      }
    }else{
      if(attributes(object)$inter == TRUE){
        selected <- which(a$all.sel[, i] != 0)
        names.i <- grep(paste0(":", a$inames[1]), rownames(a$all.sel)[selected], value = TRUE)
        names.i <- gsub(":", ".", names.i)
        if(a$boot == TRUE)
          names.i <- unique(c(names.i, grep(paste0(".", a$inames[1]), colnames(a$mat.boot[[i]]), value = TRUE)))
        if(length(names.i > 0)){
          newdata.i <- as.matrix(newdata[, gsub(paste0(".", a$inames[1]), "", names.i)]) * matrix(newdata[, a$inames[1]], nrow = nrow(newdata), ncol = length(names.i))
          colnames(newdata.i) <- names.i
          newdata <- cbind(newdata, newdata.i)
        }
      }
      if(is.null(coef(a$train.fit[[i]]))){
        new.surv$surv[, i] <- rep(a$bhz["surv"], nrow(newdata))
        new.surv$lower[, i] <- rep(a$bhz["lower"], nrow(newdata))
        new.surv$upper[, i] <- rep(a$bhz["upper"], nrow(newdata))
      }else{
        surv <- survfit(a$train.fit[[i]], newdata = newdata, se.fit = !a$boot, conf.int = a$ci.level)
        new.hrz <- max(which(surv$time < a$hrz))
        new.surv$surv[, i] <- surv$surv[new.hrz, ]
        if(a$boot == FALSE){
          new.surv$lower[, i] <- surv$lower[new.hrz, ]
          new.surv$upper[, i] <- surv$upper[new.hrz, ]
        }else{
          lpred <- a$mat.boot[[i]][, -1] %*% t(newdata[, colnames(a$mat.boot[[i]])[-1], drop = FALSE])
          newboot.surv <- exp(-matrix(a$mat.boot[[i]][, 1], nrow = nrow(lpred), ncol = ncol(lpred)) * exp(lpred))
          new.surv$lower[, i] <- apply(X = newboot.surv, MARGIN = 2, FUN = function(X) quantile(X, (1 - a$ci.level) / 2))
          new.surv$upper[, i] <- apply(X = newboot.surv, MARGIN = 2, FUN = function(X) quantile(X, 1 - (1 -a$ci.level) / 2))
        }
      }
    }
  }

  tt <- NULL
  if(a$tnames[1] == "tt")
    tt = newdata[, a$inames[1]]

  attributes(new.surv) <- append(
    x = attributes(new.surv),
    values = list(
      m.score = m.score,
      catm.score = catm.score,
      i.score = i.score,
      inter = a$inter,
      smooth = a$smooth,
      tt = tt
    ))
  class(new.surv) <- 'resexpSurv'

  return(new.surv)

} # end of predict
###########################################################

plot.resexpSurv <- function(x, method, pr.group, print.ci, xlim, ylim, xlab, ylab, ...){

  a <- attributes(x)

  if(missing(method)){
    stop("\n 'method' must be specified.")
  }else{
    if(length(method) > 1){
      warning("\n Only the first element in 'method' is considered.")
      method <- method[1]
    }
    if(!(method %in% colnames(x$surv)))
      stop("\n 'method' is unknown.")
  }

  if(a$inter == TRUE){
    if(missing(pr.group)){
      pr.group <- unique(a$catm.score[, method])
      if(a$smooth == TRUE){
        stop("\n For a meaningful visualization, 'pr.group' may be specified.")
      }else{
        warning("\n For a meaningful visualization, 'pr.group' may be specified.")
      }

    }else{
      if(sum(1 - (pr.group %in% unique(a$catm.score[, method]))) > 0)
        stop("\n 'pr.group' is not an existing group.")
    }
    wgr <- which(a$catm.score[, method] %in% pr.group)
  }else{
    if(!missing(pr.group))
      warning("\n 'pr.group' is not considered in the prognostic setting.")
    wgr <- 1:nrow(x$surv)
  }


  if(missing(print.ci)){
    print.ci <- TRUE
  }else{
    if(length(print.ci) > 1 || !is.logical(print.ci))
      stop("\n 'print.ci' must be logical.")
  }


  if(a$inter == TRUE)
    xx <- a$i.score[, method][wgr]
  else
    xx <- a$m.score[, method]

  col <- rep(1, nrow(x$surv))
  if(a$inter == TRUE)
    col <- as.numeric(factor(a$tt))[wgr]

 if(missing('xlim')) xlim <- c(0, 1)
 if(missing('ylim')) ylim <- c(0, 1)
 if(missing('ylab')) ylab <- "Survival probability"
 if(missing('xlab')) xlab <- c("Prognostic score", "Treatment-effect modifying score")[1 + a$inter]

 if(a$smooth == TRUE){
   wc <- which(col == 1)[order(xx[which(col == 1)])]
   wt <- which(col == 2)[order(xx[which(col == 2)])]
   plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
   lines(x = xx[wc], y = x$surv[wgr, method][wc], lwd = 2, col = 1)
   lines(x = xx[wt], y = x$surv[wgr, method][wt], lwd = 2, col = 2)
   if(!is.null(a$surv)){
     lines(x = xx, y = a$surv[wgr, method], type = "p", col = col, pch = "+", cex = 1)
   }
   if(print.ci){
     lines(x = xx[wc], y = x$lower[wgr, method][wc], lwd = 2, lty = 5, col = 1)
     lines(x = xx[wc], y = x$upper[wgr, method][wc], lwd = 2, lty = 5, col = 1)
     lines(x = xx[wt], y = x$lower[wgr, method][wt], lwd = 2, lty = 5, col = 2)
     lines(x = xx[wt], y = x$upper[wgr, method][wt], lwd = 2, lty = 5, col = 2)
   }
 }else{
   plot(x = xx, y = x$surv[wgr, method], col = col,
        xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
   if(print.ci){
     zz <- lapply(
       X = 1:nrow(x$surv),
       FUN = function(X){
         lines(x = rep(xx[X], 2), y = c(x$lower[wgr, method][X], x$upper[wgr, method][X]), col = col[X])
       }
     )
   }
 }
} # end of plot
###########################################################
###########################################################

int.expSSurv <- function(x, y, z, c){
  if(c == "decrease")
    pointwise <- rbind(c(-1, min(x), 1), c(1, max(x), 0))
  if(c == "increase")
    pointwise <- rbind(c(1, min(x), 0), c(-1, max(x), 1))

  sfit <- cobs(x = x, y = y, print.mesg = F, print.warn = F, method = "uniform", constraint = c, pointwise = pointwise)
  sexp <- predict(sfit, z = z)[, 2]
  sparam <-  unlist(sfit[c('knots', 'coef')])
  return(list(sexp, sparam))
}
###########################################################

int.expSurv <- function(res, all.sel, fit, datafit, dataest, boot, ci.level, hrz, out.fit){

  if(missing(fit)) nmeth <- ncol(all.sel) else nmeth <- length(fit)

  Surv <- lapply(
    X = 1:nmeth,
    FUN = function(X){

      if(out.fit){
        selected <- which(all.sel[, X] != 0)
        rownames(all.sel) <- gsub(paste0(":", attributes(res)$inames[1]), paste0(".", attributes(res)$inames[1]), rownames(all.sel))
      }

      if(attributes(res)$inter == TRUE){
        if(out.fit){
          names.i <- grep(paste0(".", attributes(res)$inames[1]), rownames(all.sel)[selected], value = TRUE)
        }else{
          names.i <- grep(paste0(".", attributes(res)$inames[1]), names(coef(fit[[X]])), value = TRUE)
        }

        if(length(names.i) > 0){
          datafit.i <- as.matrix(datafit[, gsub(paste0(".", attributes(res)$inames[1]), "", names.i)]) * matrix(datafit[, attributes(res)$inames[1]], nrow = nrow(datafit), ncol = length(names.i))
          dataest.i <- as.matrix(dataest[, gsub(paste0(".", attributes(res)$inames[1]), "", names.i)]) * matrix(dataest[, attributes(res)$inames[1]], nrow = nrow(dataest), ncol = length(names.i))
          colnames(datafit.i) <- colnames(dataest.i) <- names.i
          datafit <- cbind(datafit, datafit.i)
          dataest <- cbind(dataest, dataest.i)
        }
      }

      if(sum(selected) > 0){

        if(out.fit){

          # if(sum(selected) == 1){
            fct <- as.formula(paste0("Surv(", paste(
              attributes(res)$inames[attributes(res)$tnames == "y"], collapse = ", "),
              ") ~ ", paste0(rownames(all.sel)[selected], collapse = "+")))
          # }else{
          #   fct <- as.formula(paste0("Surv(", paste(
          #     attributes(res)$inames[attributes(res)$tnames == "y"], collapse = ", "),
          #     ") ~ as.matrix(datafit[, rownames(all.sel)[selected], drop = FALSE])"))
          # }

          fit <- coxph(fct, init = all.sel[selected, X], iter = 0, data = datafit)
        }
          surv <- survfit(fit, newdata = dataest, se.fit = !boot, conf.int = ci.level)
          whrz <- max(which(surv$time < hrz))
          if(boot == TRUE){
            surv <- list(surv = surv$surv[whrz, ])
          }else{
            surv <- list(surv = surv$surv[whrz, ], lower = surv$lower[whrz, ], upper = surv$upper[whrz, ])
          }
      }else{
        fct <- as.formula(paste0("Surv(", paste(
          attributes(res)$inames[attributes(res)$tnames == "y"], collapse = ", "),
          ") ~ 1"))
        fit <- coxph(fct, data = dataest)
        surv <- survfit(fct, data = dataest, se.fit = !boot, conf.int = ci.level)
        whrz <- max(which(surv$time < hrz))
        if(boot == TRUE){
          surv <- list(surv = rep(surv$surv[whrz], length(surv$surv)))
        }else{
          surv <- list(surv = rep(surv$surv[whrz], length(surv$surv)),
                       lower = rep(surv$lower[whrz], length(surv$surv)),
                       upper = rep(surv$upper[whrz], length(surv$surv)))
        }
      }

      if(out.fit){
        return(list(surv, fit))
      }else{
        return(list(surv))
      }
    })
  return(Surv)
}
###########################################################

int.predSmooth <- function(i, label, score, wrow, mat){
  ocobs <- list(tau = .5, degree = 2)
  class(ocobs) <- "cobs"
  if(missing(wrow))
    wrow = 1:nrow(score)
  ocobs$knots <- setdiff(mat[label, 1:6], 0)
  ocobs$coef <- setdiff(mat[label, 7:13], 0)
  predict(ocobs, z = score[wrow, i])[, 2]
}
###########################################################
