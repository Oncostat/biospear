################################################################################
### Elastic-net penalty                                                        #
################################################################################

enet <- function(data, o){
  
  x <- as.matrix(data[, attributes(data)$pos$X])
  y <- as.matrix(data[, attributes(data)$pos$y])
  
  ### Cross-validation to estimate the optimal (alpha, lambda) combination
  cv <- matrix(unlist(lapply(
    X = seq(0.05, 0.95, 0.05),
    FUN = function(X){
      
      lGrid <- l.Grid(
        x = x,
        y = y,
        w = attributes(data)$weights,
        nfold = length(unique(o$folds)),
        alpha = X,
        dfmax = o$dfmax,
        family = "cox")
      
      cv <- cv.glmnet(
        x = x,
        y = y,
        family = "cox",
        alpha = X,
        lambda = seq(lGrid[1], lGrid[2], diff(lGrid) / 100),
        penalty.factor = attributes(data)$weights,
        foldid = o$folds,
        grouped = TRUE,
        standardize = FALSE,
        thresh = 1e-16)
      res <- c(cv$cvm[which(cv$lambda == cv$lambda.min)], cv$lambda.min, X)
      return(res)
    })), ncol = 3, byrow = T, dimnames = list(NULL, c("cvm", "lambda", "alpha")))
  
  ### Sequential order of removing biomarkers
  step <- NA
  if(o$isSim == TRUE){
    step <- rank(rowSums(
      glmnet(
        x = x,
        y = y,
        family = "cox",
        alpha = cv[which.max(cv[, "cvm"]), "alpha"],
        penalty.factor = attributes(data)$weights,
        nlambda = ncol(x),
        standardize = FALSE,
        thresh = 1e-16)$beta != 0
    )[which(attributes(data)$weights != 0)], ties.method = "min")
  }
  
  ### Fit of the final model
  nfold <- length(unique(o$folds))
  fit.enet <- glmnet(
    x = x,
    y = y,
    family = "cox",
    penalty.factor = attributes(data)$weights,
    alpha = cv[which.max(cv[, "cvm"]), "alpha"],
    lambda = cv[which.max(cv[, "cvm"]), "lambda"] * ((nfold - 1) / nfold),
    standardize = FALSE,
    thresh = 1e-16)
  
  names <- coef(fit.enet)@Dimnames[[1]][which(coef(fit.enet) != 0)]
  res <- coef(fit.enet)[which(coef(fit.enet) != 0)]
  names(res) <- names
  
  return(list(res, step))
}