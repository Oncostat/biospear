################################################################################
### PCA + lasso (only available for predictive setting)                        #
################################################################################

PCAlasso <- function(data, o){

  x <- as.matrix(data[, attributes(data)$pos$x])
  x_ <- as.matrix(data[, unlist(attributes(data)$pos[c('tt', 'z', 'xt')])])
  xt <- as.matrix(data[, attributes(data)$pos$xt])
  y <- as.matrix(data[, attributes(data)$pos$y])
  
  #########################################################
  ### Preliminary step (dimension reduction of main effects matrix via PCA)
  
  ## Estimation of the ncompS first principal components
  ncompS <- round(min(ncol(x),
                      table(data[which(data$status == 1), "treat"])))

  
  pcS <- (fast.svd(t(t(x) - colMeans(x))))
  PCS <- as.matrix(
    pcS$u[, 1:ncompS, drop = FALSE] %*% diag(pcS$d[1:ncompS], ncompS))
  colnames(PCS) <- gsub(" ", "", format(paste0("PCS", 1:ncompS)))
  beta <- pcS$v
  PCfolds <- sample(o$folds) 
  
  ## Cross-validation to estimate the optimal number K
  
  # if(o$family != "cox"){
  #   pos.y <- grep(c("y.cont", "y.bin")[which(c("gaussian", "binomial") %in% o$family)], names(data))
  #   cv.K <- t(sapply(1:max(PCfolds),FUN = function(i) {
  #     PCS.T <- PCS[which(PCfolds != i),]
  #     sapply(X = 1:ncompS, FUN = function(X){  
  #       fit <- glm(data[, pos.y] ~ PCS.T[,1:X],
  #                    data = data[which(PCfolds != i),])
  #       coef.T <- fit$coefficients
  #       nlogl.T <- fit$loglik[2]
  #       off.All <- PCS[, 1:X] %*% as.matrix(coef.T) 
  #       nlogl.All <- glm(data[, pos.y] ~ offset(off.All),
  #                          data = data)$loglik
  #       nlogl <- nlogl.All - nlogl.T
  #       return(nlogl)
  #     })
  #   }))    
  # }
  
  cv.K <- t(sapply(1:max(PCfolds),FUN = function(i) {
    PCS.T <- PCS[which(PCfolds != i), ]
    sapply(X = 1:ncompS, FUN = function(X){  
      fit <- coxph(Surv(time, status) ~ PCS.T[,1:X],
                   data = data[which(PCfolds != i), ])
      coef.T <- fit$coefficients
      nlogl.T <- fit$loglik[2]
      off.All <- PCS[, 1:X] %*% as.matrix(coef.T) 
      nlogl.All <- coxph(Surv(time, status) ~ offset(off.All),
                         data = data)$loglik
      nlogl <- nlogl.All - nlogl.T
      return(nlogl)
    })
  }))

  K.PCA <- which.max(colSums(cv.K))
  x_ <- as.matrix(cbind(x_, data.frame(PCS[, 1:K.PCA])))
  beta <- beta[, 1:K.PCA]
  pos.PCS <- grep("PCS", names(x_))

  #########################################################
  ### Final model
  
  w <- rep(0, ncol(x_))
  w[c(grep("bm", colnames(x_)), grep("bi", colnames(x_)))] <- 1
  
  nfold <- length(unique(o$folds))
  lGrid <- l.Grid(
    x = x_,
    y = y,
    w = w,
    nfold = nfold,
    dfmax = o$dfmax,
    family = "cox")
  
  ### Cross-validation to estimate the optimal lambda
  cv <- cv.glmnet(
    x = x_,
    y = y,
    family = "cox",
    alpha = 1,
    lambda = seq(lGrid[1], lGrid[2], diff(lGrid) / 100),
    penalty.factor = w,
    foldid = o$folds,
    standardize = FALSE,
    thresh = 1e-16)
  l <- cv$lambda.min
  
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
    lambda = l * ((nfold - 1) / nfold),
    penalty.factor = w,
    standardize = FALSE,
    thresh = 1e-16)

  names <- coef(fit)@Dimnames[[1]][which(coef(fit) != 0)]
  res <- coef(fit)[which(coef(fit) != 0)]
  names(res) <- names
  res <- res[-grep("PCS", names(res))]
  res.x <- as.vector(beta %*% t(t(coef(fit)[grep("PCS", rownames(coef(fit)))])))
  names(res.x) <- colnames(x)
  res <- c(res, res.x)
  
  return(list(res, step))
}