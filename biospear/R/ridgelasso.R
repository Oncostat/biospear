################################################################################
### Offset Ridge+lasso penalty (only available for predictive setting)         #
################################################################################

ridgelasso <- function(data, o){

  x <- as.matrix(data[, attributes(data)$pos$x])
  x_ <- as.matrix(data[, unlist(attributes(data)$pos[c('tt', 'z', 'xt')])])
  xt <- as.matrix(data[, attributes(data)$pos$xt])
  y <- as.matrix(data[, attributes(data)$pos$y])
  
  #########################################################
  ### Preliminary step to estimate off.x
  
  ## Cross-validation to estimate the optimal lambda
  folds <- sample(o$folds)
  cv.x <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 0,
    foldid = folds,
    grouped = TRUE,
    standardize = FALSE,
    thresh = 1e-16)
  
  ## Fit of the model
  nfold <- length(unique(o$folds))
  fit.x <- glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 0,
    lambda = cv.x$lambda.min * ((nfold - 1) / nfold),
    standardize = FALSE,
    thresh = 1e-16)
  res.x <- as.vector(coef(fit.x))
  
  ## Offset for the final model
  off.x <- x %*% as.matrix(res.x)
  names(res.x) <- rownames(coef(fit.x))
  
  #########################################################
  
  #########################################################
  ### Final model
  
  ## Cross-validation to estimate the optimal lambda
  lGrid <- l.Grid(
    x = x_,
    y = y,
    w = attributes(data)$weights[-c(length(attributes(data)$weights) - (0:(ncol(x) - 1)))],
    nfold = length(unique(o$folds)),
    off = off.x,
    dfmax = o$dfmax,
    family = "cox")
  
  cv.x_ <- cv.glmnet(
    x = x_,
    y = y,
    family = "cox",
    alpha = 1,
    lambda = seq(lGrid[1], lGrid[2], diff(lGrid) / 100),
    foldid = o$folds,
    grouped = TRUE,
    offset = off.x,
    penalty.factor = attributes(data)$weights[-c(length(attributes(data)$weights) - (0:(ncol(x) - 1)))],
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
        offset = off.x,
        penalty.factor = attributes(data)$weights[-c(length(attributes(data)$weights) - (0:(ncol(x) - 1)))],
        nlambda = ncol(x),
        standardize = FALSE,
        thresh = 1e-16)$beta != 0
    )[which(attributes(data)$weights[-c(length(attributes(data)$weights) - (0:(ncol(x) - 1)))] != 0)], ties.method = "min")
    step.x <- rep(NA, ncol(x))
    names(step.x) <- colnames(x)
    step <- c(step.x, step)
  }
  
  ## Fit of the model
  fit.x_ <- glmnet(
    x = x_,
    y = y,
    family = "cox",
    alpha = 1,
    lambda = cv.x_$lambda.min * ((nfold - 1) / nfold),
    offset = off.x,
    penalty.factor = attributes(data)$weights[-c(length(attributes(data)$weights) - (0:(ncol(x) - 1)))],
    standardize = FALSE,
    thresh = 1e-16)
  
  names <- coef(fit.x_)@Dimnames[[1]][which(coef(fit.x_) != 0)]
  res <- coef(fit.x_)[which(coef(fit.x_) != 0)]
  names(res) <- names
  res <- c(res, res.x)
  
  return(list(res, step))
}