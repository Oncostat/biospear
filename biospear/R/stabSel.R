################################################################################
### Stability selection                                                        #
################################################################################

stabSel <- function(data, o){
  
  x <- as.matrix(data[, attributes(data)$pos$X])
  y <- as.matrix(data[, attributes(data)$pos$y])
  
  fit <- glmnet(
    x = x,
    y = y,
    family = "cox",
    penalty.factor = attributes(data)$weights,
    alpha = 1,
    nlambda = ncol(x))
  
  subsets <- sapply(1:o$ss.nsub, function(v) {
    sample(1:nrow(x), nrow(x) * o$ss.fsub)
  })
  
  res.ss <- lapply(
    X = 1:o$ss.nsub,
    FUN = function(X){
      w <- attributes(data)$weights
      if(o$ss.rando == TRUE)
        w[which(w == 1)] <- 1/runif(length(which(w == 1)), 0.5, 1)
      fit.ss <- glmnet(x = x[subsets[, X], ],
                    y = y[subsets[, X], ],
                    family = "cox",
                    penalty.factor = w,
                    alpha = 1,
                    lambda = fit$lambda)
      res <- as.matrix(coef(fit.ss))
      res <- (res != 0) * 1
      return(res)
    })
  freq.biom <- matrix(0, nrow = nrow(res.ss[[1]]), ncol = ncol(res.ss[[1]]))
  for(i in seq(length(res.ss)))
    freq.biom <-  freq.biom + (res.ss[[i]]/o$ss.nsub)
  if(o$inter == TRUE){
    which.sel <- grep("bi", rownames(freq.biom))    
  }else{
    which.sel <- grep("bm", rownames(freq.biom))    
  }
  nb.biom <- colSums(freq.biom[which.sel, ])
  
  qv <- ceiling(sqrt(o$ss.fwer * (2 * o$ss.thr - 1) * length(which.sel)))

  if(nb.biom[length(nb.biom)] <= qv){
    lpos <- length(length(nb.biom))
  }
  else{
    lpos <- which(nb.biom > qv)[1]
  }
  
  stable <- which(freq.biom[, lpos] >= o$ss.thr)
  ffit <- as.formula(paste0("Surv(time, status) ~ ", paste0(names(stable), collapse = "+")))
  res <- coef(coxph(ffit, data = data))
  
  step <- rowSums(freq.biom >= o$ss.thr)[which(attributes(data)$weights != 0)]
  
  return(list(res, step))
}