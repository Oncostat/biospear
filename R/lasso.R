################################################################################
### Lasso penalty                                                              #
################################################################################

lasso <- function(data, o){
  
  x <- as.matrix(data[, attributes(data)$pos$X])
  y <- as.matrix(data[, attributes(data)$pos$y])
  
  lGrid <- l.Grid(
    x = x,
    y = y,
    w = attributes(data)$weights,
    nfold = length(unique(o$folds)),
    dfmax = o$dfmax,
    family = "cox")
  
  ### Cross-validation to estimate the optimal lambda
  cv <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    lambda = seq(lGrid[1], lGrid[2], diff(lGrid) / 100),
    foldid = o$folds,
    grouped = TRUE,
    penalty.factor = attributes(data)$weights,
    standardize = FALSE,
    thresh = 1e-16)
  l <- cv$lambda.min
  
  ### Sequential order of removing biomarkers
  step <- NA
  if(o$isSim == TRUE){
    step <- rank(rowSums(
      glmnet(
        x = x,
        y = y,
        family = "cox",
        alpha = 1,
        penalty.factor = attributes(data)$weights,
        nlambda = ncol(x),
        standardize = FALSE,
        thresh = 1e-16)$beta != 0
    )[which(attributes(data)$weights != 0)], ties.method = "min")
  }
  
  ### Fit of the models for different lambda values
  nfold <- length(unique(o$folds))
  fit.lasso <- glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    penalty.factor = attributes(data)$weights,
    lambda = l * ((nfold - 1) / nfold),
    standardize = FALSE,
    thresh = 1e-16)
  
  names <- coef(fit.lasso)@Dimnames[[1]][which(coef(fit.lasso) != 0)]
  res <- coef(fit.lasso)[which(coef(fit.lasso) != 0)]
  names(res) <- names

  ### Save for Adaptive lasso (preliminary step with lasso, alassoL)
  pr.alsL <- coef(fit.lasso)[-c(1:(sum(attributes(data)$weights == 0)))]
  
  ### Save for empirical penalized lasso extensions (penExt)
    fit <- glmnet(
      x = x,
      y = y,
      family = "cox",
      alpha = 1,
      penalty.factor = attributes(data)$weights,
      lambda = cv$lambda * ((nfold - 1) / nfold),
      standardize = FALSE,
      thresh = 1e-16)
  
  return(list(res, 
              step,
              l,
              pr.alsL,
              cv,
              fit
              )
         )
}