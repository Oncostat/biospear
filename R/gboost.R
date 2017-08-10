################################################################################
### Gradient boosting                                                          #
################################################################################

gboost <- function(data, o){

  x <- as.matrix(data[, attributes(data)$pos$X])
  y <- as.matrix(data[, attributes(data)$pos$y])
  
  ### Treatment and clinical covariates effects (offset)
  off.zt <- res.zt <- NULL
  if(sum(attributes(data)$weights == 0) > 0){
    
    f <- as.formula(paste("Surv(time, status) ~ ", 
                          paste(names(data)[c(attributes(data)$pos$t, attributes(data)$pos$z)], collapse = "+")))
    fit.zt <- coxph(f, data = data)
    res.zt <- coef(fit.zt)
    off.zt <- fit.zt$linear.predictors

    
    # if(o$family != "cox"){
    #   yf <- c("y.cont", "y.bin")[which(c("gaussian", "binomial") %in% o$family)]
    #   
    #   f <- as.formula(paste(yf, 
    #                         paste(names(data)[c(attributes(data)$pos$t, attributes(data)$pos$z)], collapse = "+")))
    #   fit.zt <- glm(f, family = o$family, data = data)
    #   res.zt <- coef(fit.zt)
    #   off.zt <- fit.zt$linear.predictors
    # } 
      
  }
  
  # Gradient boosting
  y <- Surv(y[,1], y[,2])
  fam <- CoxPH()
 
  # if(o$family == "gaussian") fam <- Gaussian()
  # if(o$family == "binomial") fam <- Binomial() 
  
  boost <- glmboost(
      x = x[, c(grep("bm", colnames(x)), grep("bi", colnames(x)))],
      y = y,
      offset = off.zt,
      family = fam,
      center = FALSE,
      control = boost_control(mstop = 2 * ncol(x)))

  step <- NA
  if(o$isSim == TRUE){
    namesbiom <- unlist(lapply(
      X = 1:1000,
      FUN = function(X)
        names(coef(boost[X]))
    ))
    step <- table(
      c(namesbiom, 
        colnames(data[, c(attributes(data)$pos$x, attributes(data)$pos$xt)]))
    ) - 1
    step <- c(step[-(1:(length(step)/2))], step[(1:(length(step)/2))])
  }
  
  # Cross-validation to find the optimal number of steps
  nfold <- length(unique(o$folds))
  cv <- cv(
    weights = model.weights(boost),
    type = "kfold", 
    B = nfold)
  
  cv <- cvrisk(
    object = boost,
    grid = 1:mstop(boost),
    folds = cv,
    papply = lapply)
  
  # Fit of the final model
  fit <- boost[mstop(cv)]
  
  res <- coef(fit)
  attributes(res) <- NULL
  names(res) <- names(coef(fit))
  res <- c(res.zt, res)
  
  return(list(res, step))
}
