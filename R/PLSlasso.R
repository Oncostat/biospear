################################################################################
### PLS + lasso (only available for predictive setting)                        #
################################################################################

PLSlasso <- function(data, o){

  x <- as.matrix(data[, attributes(data)$pos$x])
  x_ <- as.matrix(data[, unlist(attributes(data)$pos[c('tt', 'z', 'xt')])])
  xt <- as.matrix(data[, attributes(data)$pos$xt])
  y <- as.matrix(data[, attributes(data)$pos$y])

  #########################################################
  ### Preliminary step (dimension reduction of main effects matrix via PLS)
  ### Estimation of the first component

  pls <- coxpls(
    Xplan = x,
    time = as.matrix(y[, 1]),
    event = as.matrix(y[, 2]),
    ncomp = 1,
    allres = T)
  PLS <- as.matrix(pls$tt_pls)
  x_ <- cbind(x_, PLS)
  pos.PLS <- grep("dim.1", colnames(x_))

  #########################################################
  ### Final model

  ### Lambda.max estimation
  w <- rep(0, ncol(x_))
  w[grep("bi", colnames(x_))] <- 1

  nfold <- length(unique(o$folds))
  lGrid <- l.Grid(
    x = x_,
    y = y,
    w = w,
    nfold = nfold,
    dfmax = o$dfmax,
    family = "cox")

  ### Cross-validation to estimate the optimal lambda
  cv.PLS <- cv.glmnet(
    x = x_,
    y = y,
    family = "cox",
    alpha = 1,
    lambda = seq(lGrid[1], lGrid[2], diff(lGrid) / 100),
    penalty.factor = w,
    foldid = o$folds,
    standardize = FALSE,
    thresh = 1e-16)

  ### Sequential order of removing biomarkers
  step <- NA
  if(o$isSim == TRUE){
    step <- rank(rowSums(
      glmnet(
        x = x_,
        y = y,
        family = "cox",
        alpha = 1,
        penalty.factor = w,
        nlambda = ncol(x) * 2,
        standardize = FALSE,
        thresh = 1e-16)$beta != 0
    )[which(w != 0)], ties.method = "min")
    step.x <- rep(NA, ncol(x))
    names(step.x) <- colnames(x)
    step <- c(step.x, step)
  }

  ### Fit of the final model
  fit <- glmnet(
    x = x_,
    y = y,
    family = "cox",
    alpha = 1,
    lambda = cv.PLS$lambda.min * ((nfold - 1) / nfold),
    penalty.factor = w,
    standardize = FALSE,
    thresh = 1e-16)

  names <- coef(fit)@Dimnames[[1]][which(coef(fit) != 0)]
  res <- coef(fit)[which(coef(fit) != 0)]
  names(res) <- names
  res <- res[-length(res)]

  res.x <- as.vector(as.matrix(pls$pls_mod$loadings$X) %*% t(t(coef(fit)[nrow(coef(fit))])))
  names(res.x) <- colnames(x)
  res <- c(res, res.x)

  return(list(res, step))
}
