################################################################################
### Modified covariates (Tian et al., JASA 2014)                               #
################################################################################

modCov <- function(data, o){
  
  x <- as.matrix(data[, attributes(data)$pos$x])
  z <- as.matrix(data[, attributes(data)$pos$z])
  xt <- as.matrix(scale(data[, attributes(data)$pos$xt], center = T, scale = T))
  zxt <- as.matrix(cbind(z, xt))
  y <- as.matrix(data[, attributes(data)$pos$y])
  w <- c(rep(0, ncol(z)), rep(1, ncol(xt)))
  
  #########################################################
  ### Final model
  
  ## Cross-validation to estimate the optimal lambda
  nfold <- length(unique(o$folds))
  lGrid <- l.Grid(
    x = zxt,
    y = y,
    w = w,
    nfold = nfold,
    dfmax = o$dfmax,
    family = "cox")
  
  cv.modCov <- cv.glmnet(
    x = zxt,
    y = y,
    family = "cox",
    alpha = 1,
    lambda = seq(lGrid[1], lGrid[2], diff(lGrid) / 100),
    foldid = o$folds,
    grouped = TRUE,
    penalty.factor = w,
    standardize = FALSE,
    thresh = 1e-16)
  
  ### Sequential order of removing biomarkers
  step <- NA
  if(o$isSim == TRUE){
    step <- rank(rowSums(
      glmnet(
        x = zxt,
        y = y,
        family = "cox",
        alpha = 1,
        penalty.factor = w,
        nlambda = ncol(xt) * 2,
        standardize = FALSE,
        thresh = 1e-16)$beta != 0
    )[which(w != 0)], ties.method = "min")
    step.x <- rep(NA, ncol(x))
    names(step.x) <- colnames(x)
    step <- c(step.x, step)
  }
  
  ## Fit of the model
  fit.modCov <- glmnet(
    x = zxt,
    y = y,
    family = "cox",
    alpha = 1,
    lambda = cv.modCov$lambda.min * ((nfold - 1) / nfold),
    penalty.factor = w,
    standardize = FALSE,
    thresh = 1e-16)
  
  names <- coef(fit.modCov)@Dimnames[[1]][which(coef(fit.modCov) != 0)]
  res <- coef(fit.modCov)[which(coef(fit.modCov) != 0)]
  names(res) <- names
  
  return(list(res, step))
}