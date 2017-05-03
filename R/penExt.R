################################################################################
### Empirical penalized extensions of the lasso penalty                        #
################################################################################

penExt <- function(data, o){

  x <- as.matrix(data[, attributes(data)$pos$X])
  y <- as.matrix(data[, attributes(data)$pos$y])  
  nfold <- length(unique(o$folds))
  
  cv <- unlist(o$cv.lasso, recursive = FALSE)
  fit <- unlist(o$fit.lasso, recursive = FALSE)
  l.cvl <- unlist(o$l.lasso, recursive = FALSE)

  # Non-zero parameters (nz) and cross-validated loglikelihood (cvl) for all the computed lambdas
  if(o$inter == FALSE){
    nz <- rev(fit$df) - sum(attributes(data)$weight == 0)
  }else{
    nz <- rev(colSums(fit$beta[attributes(data)$pos$xt - 2, ] != 0))
  }
  l <- rev(cv$glmnet.fit$lambda[1 : length(fit$lambda)])
  cvl <- rev((cv$cvm*(-nrow(data)/2))[1:length(fit$lambda)])
  cvl0 <- max(cvl[which(nz == 0)])
  cvlMax <- cvl[which(l == l.cvl)]
  nzMax <- nz[which(l == l.cvl)]

  l.Ext <- rep(l.cvl, length(o$n.penExt))
  
  if(nzMax != 0){
    
    # lasso-pcvl: penalized cross-validated loglikelihood (Ternes et al., 2016)
    l.pcvl <- NULL
    if("lasso-pcvl" %in% o$n.penExt){
      pcvl <- cvl - (((cvlMax - cvl0) / nzMax) * nz)
      l.pcvl <- l[which.max(pcvl)]
      names(l.pcvl) <- "lasso-pcvl"
    }
    
    # lasso-AIC: Akaike Information Criterion
    l.AIC <- NULL
    if("lasso-AIC" %in% o$n.penExt){
      aic <- cvl - (2 * nz)
      l.AIC <- l[which.max(aic)]
      names(l.AIC) <- "lasso-AIC"
    }
    
    # lasso-BIC: Bayesian Information Criterion
    l.BIC <- NULL
    if("lasso-BIC" %in% o$n.penExt){
      bic <- cvl - (log(nrow(data)) * nz)
      l.BIC <- l[which.max(bic)]
      names(l.BIC) <- "lasso-BIC"
    }
    
    # lasso-AICC: corrected AIC
    l.AICC <- NULL
    if("lasso-AICC" %in% o$n.penExt){
      aicc <- cvl - (2 * ((nz + 1) / (nrow(data) - nz - 2)))
      l.AICC <- l[which.max(aicc)]
      names(l.AICC) <- "lasso-AICC"
    }
    
    # lasso-RIC: Risk Information Criterion
    l.RIC <- NULL
    if("lasso-RIC" %in% o$n.penExt){
      ric <- cvl - (2 * log(nz) * nz)
      l.RIC <- l[which.max(ric)]
      names(l.RIC) <- "lasso-RIC"
    }
    
    # lasso-HQIC: Hannan and Quinn Information Criterion
    l.HQIC <- NULL
    if("lasso-HQIC" %in% o$n.penExt){
      hqic <- cvl - (2 * nz * log(log(nrow(data))))
      l.HQIC <- l[which.max(hqic)]
      names(l.HQIC) <- "lasso-HQIC"
    }
    
    # lasso-eBIC: extented BIC
    l.eBIC <- NULL
    if("lasso-eBIC" %in% o$n.penExt){
      ebic <- cvl - ((log(nrow(data)) * nz) + (2 * 0.5 * nz * log(ncol(x))))
      l.eBIC <- l[which.max(ebic)]
      names(l.eBIC) <- "lasso-eBIC"
    }

    # lasso-1se : one-standard-error rule (Friedman et al., 2010)
    l.1se <- NULL
    if("lasso-1se" %in% o$n.penExt){
      l.1se <- cv$lambda.1se
      names(l.1se) <- "lasso-1se"
    }
    
    # lasso-pct : percentile lasso (Roberts and Nowak, 2013)
    l.pct <- NULL
    if("lasso-pct" %in% o$n.penExt){
      l.pct <- quantile(unlist(lapply(
        X = 1:o$pct.rep,
        FUN = function(X){
          folds <- sample(o$folds)
          lGrid <- l.Grid(
            x = x, 
            y = y,
            w = attributes(data)$weights,
            family = "cox",
            nfold = length(unique(folds)),
            dfmax = o$dfmax)
      
          ### Cross-validation to estimate the optimal lambda
          l <- cv.glmnet(
            x = x,
            y = y,
            family = "cox",
            alpha = 1,
            lambda = seq(lGrid[1], lGrid[2], diff(lGrid) / 100),
            foldid = folds,
            penalty.factor = attributes(data)$weights,
            grouped = TRUE,
            standardize = FALSE,
            thresh = 1e-16)$lambda.min
          return(l)
        })), probs = o$pct.qtl)
      names(l.pct) <- "lasso-pct"
    }
    l.Ext <- c(l.pcvl, l.AIC, l.BIC, l.AICC, l.RIC, l.HQIC, l.eBIC, l.1se, l.pct)
  }

    res <- lapply(
    X = 1:length(l.Ext),
    FUN = function(X){
      fit.Ext <- glmnet(
        x = x,
        y = y,
        family = "cox",
        alpha = 1,
        lambda = l.Ext[X] * ((nfold - 1) / nfold),
        penalty.factor = attributes(data)$weights,
        standardize = FALSE,
        thresh = 1e-16)
      
      names.Ext <- coef(fit.Ext)@Dimnames[[1]][which(coef(fit.Ext) != 0)]
      res.Ext <- coef(fit.Ext)[which(coef(fit.Ext) != 0)]
      names(res.Ext) <- names.Ext
      return(res.Ext)
    })

  names(res) <- names(l.Ext)
  
  return(res)
}