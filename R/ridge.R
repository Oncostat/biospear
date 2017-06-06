################################################################################
### Ridge penalty                                                              #
################################################################################

ridge <- function(data, o){
  
  x <- as.matrix(data[, attributes(data)$pos$X])
  y <- as.matrix(data[, attributes(data)$pos$y])
  
  ### Cross-validation to estimate the optimal lambda
  cv <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 0,
    foldid = o$folds,
    penalty.factor = attributes(data)$weights,
    grouped = TRUE,
    standardize = FALSE,
    thresh = 1e-16)
  l <- cv$lambda.min
  
  ### Sequential order of removing biomarkers
  step <- NA
  if(o$isSim == TRUE){
    step <- rep(0, length(unlist(attributes(data)$inames[c('x', 'xt')])))
  }
  
  ### Fit of the models for different lambda values
  nfold <- length(unique(o$folds))
  fit <- glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 0,
    lambda = l * ((nfold - 1) / nfold),
    penalty.factor = attributes(data)$weights,
    standardize = FALSE,
    thresh = 1e-16)
  
  names <- coef(fit)@Dimnames[[1]][which(coef(fit) != 0)]
  res <- coef(fit)[which(coef(fit) != 0)]
  names(res) <- names

  ### Save for Adaptive lasso (preliminary step with ridge, alassoR)
  pr.alsR <- coef(fit)[-c(1:(sum(attributes(data)$weights == 0)))]
  
  return(list(res, step, pr.alsR))
}