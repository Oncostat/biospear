################################################################################
### Adaptive lasso penalty (Zou, 2006; Zhang and Lu, 2007)                     #
################################################################################

alasso <- function(data, o, pr){
  
  coef <- unlist(eval(parse(text = paste0('o$pr.als', pr))))
  attr <- attributes(data)
  if(length(coef[which(coef != 0)]) > sum(attributes(data)$weights == 0)){
    
    x <- as.matrix(data[, c(c(attr$pos$tt, attr$pos$z), 
                            c(attr$pos$x, attr$pos$xt)[which(coef != 0)])])
    y <- as.matrix(data[, attributes(data)$pos$y])
    
    coef <- coef[which(coef != 0)]
    w <- 1 / abs(coef)
    w <- w / mean(w)
    w <- c(rep(0, sum(attributes(data)$weights == 0)), w)
    
    lGrid <- l.Grid(
      x = x,
      y = y,
      w = w,
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
      standardize = FALSE,
      penalty.factor = w,
      thresh = 1e-16)
    l <- cv$lambda.min
    
    ### Sequential order of removing biomarkers
    step <- NA
    if(o$isSim == TRUE){
      step <- rep(0, sum(attributes(data)$weights != 0))
      Step <- rowSums(
        glmnet(
          x = x,
          y = y,
          family = "cox",
          alpha = 1,
          penalty.factor = w,
          nlambda = ncol(x),
          standardize = FALSE,
          thresh = 1e-16)$beta != 0
      )[which(w != 0)]
      step[which(attributes(data)$names[unlist(attributes(data)$pos[c('x', 'xt')])] %in% names(Step))] <- Step
      step <- rank(step, ties.method = "min")
    }
    
    ### Fit of the models for different lambda values
    nfold <- length(unique(o$folds))
    fit <- glmnet(
      x = x,
      y = y,
      family = "cox",
      alpha = 1,
      lambda = l * ((nfold - 1) / nfold),
      penalty.factor = w,
      standardize = FALSE,
      thresh = 1e-16)
     
  }else{
    nfold <- length(unique(o$folds))
    x <- as.matrix(data[, attributes(data)$pos$X])
    y <- as.matrix(data[, attributes(data)$pos$y])
    
    fit <- glmnet(
      x = x,
      y = y,
      family = "cox",
      alpha = 1,
      penalty.factor = attributes(data)$weights,
      lambda = o$l.lasso * ((nfold - 1) / nfold),
      standardize = FALSE,
      thresh = 1e-16)
    
  }
      
  names <- coef(fit)@Dimnames[[1]][which(coef(fit) != 0)]
  res <- coef(fit)[which(coef(fit) != 0)]
  names(res) <- names
  
  return(list(res, step))
}