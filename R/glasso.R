################################################################################
### group-lasso penalty (only available for predictive setting)                #
################################################################################

glasso <- function(data, o){

  x <- as.matrix(data[, attributes(data)$pos$X])
  y <- as.matrix(data[, attributes(data)$pos$y])
  nfold <- length(unique(o$folds))
  x <- cbind(id = 1:nrow(x), x, foldid = o$folds)

  ### Data transformation
  poidata <- poissonize(
    data = as.data.frame(cbind(x, y)),
    interval.width = 1/6,
    factors = colnames(x)[-1])
  group <- c(rep(NA, sum(attributes(data)$weight == 0)), rep(1:(sum(attributes(data)$weight == 1)/2), 2))

  ### Lambda.max estimation
  lmax <- ceiling(
    lambdamax(
      x = as.matrix(poidata[, 1 + (1:length(group))]),
      y = poidata$m,
      index = group,
      offset = log(poidata$Rt),
      center = FALSE,
      standardize = FALSE,
      data = poidata,
      model = PoissReg(),
      control = grpl.control(trace = 0)) + 2)
  lseq <- rev(seq(2, lmax, 1))

  ### Cross-validation to estimate the optimal lambda
  cv.gLasso <- t(sapply(1:nfold, FUN = function(i){
    poidata.T <- poidata[which(poidata$foldid != i), ]
    poidata.V <- poidata[which(poidata$foldid == i), ]
    res.T <- grplasso(
      x = as.matrix(poidata.T[, 1 + (1:length(group))]),
      y = poidata.T$m,
      index = group,
      offset = log(poidata.T$Rt),
      center = FALSE,
      standardize = FALSE,
      lambda = rev(lseq),
      model = PoissReg(),
      control = grpl.control(trace = 0))
    off.V <- as.matrix(poidata.V[, 1 + (1:length(group))]) %*% res.T$coefficients
    res.V <- sapply(
      X = 1:length(lseq),
      FUN = function(X){
        grplasso(
          x = as.matrix(poidata.V[, 3]),
          y = poidata.V$m,
          index = 1,
          offset = log(poidata.V$Rt) + as.vector(off.V[,X]),
          lambda = c(1e10, 0),
          center = FALSE,
          standardize = FALSE,
          model = PoissReg(),
          data = poidata.V,
          control = grpl.control(trace = 0))$nloglik[1]
        })
    return(res.V)
  }))
  cvpos.gLasso <- lseq[which.min(apply(cv.gLasso, 2, sum))]

  ### Fit of multiple models
  step <- NA
  if(o$isSim == TRUE){
    step <- grplasso(
      x = as.matrix(poidata[, 1 + (1:length(group))]),
      y = poidata$m,
      index = group,
      offset = log(poidata$Rt),
      lambda = lseq,
      center = FALSE,
      standardize = FALSE,
      data = poidata,
      model = PoissReg(),
      control = grpl.control(trace = 0))
    step <- rowSums(step$coefficients != 0)[which(attributes(data)$weight != 0)]
  }

  ### Fit of the final model
  fit.glasso <- grplasso(
    x = as.matrix(poidata[, 1 + (1:length(group))]),
    y = poidata$m,
    index = group,
    offset = log(poidata$Rt),
    lambda = cvpos.gLasso,
    center = FALSE,
    standardize = FALSE,
    data = poidata,
    model = PoissReg(),
    control = grpl.control(trace = 0))
  names <- rownames(coef(fit.glasso))[which(coef(fit.glasso) != 0)]
  res <- coef(fit.glasso)[which(coef(fit.glasso) != 0)]
  names(res) <- names

  return(list(res, step))

}
