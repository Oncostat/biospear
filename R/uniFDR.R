################################################################################
### Univariate selection approach with FDR control                             #
################################################################################

uniFDR <- function(data, o){
  
  if(o$inter == FALSE) id <- c(rep(0, sum(attributes(data)$weights == 0)), 1:sum(attributes(data)$weights == 1))    
  if(o$inter == TRUE) id <- c(rep(0, sum(attributes(data)$weights == 0)), rep(1:(sum(attributes(data)$weights == 1) / 2), 2))   
    
  # Computation of adjusted p-values (Benjamini and Hochberg, 1995)
  fwer <- p.adjust(p = o$p.uniFDR, method = "BH")
  sign <- fwer <= o$uni.fdr
  
  # Sequential order of removing biomarkers
  step <- NA
  if(o$isSim == TRUE)
    step <- rank(-(rep(fwer, ifelse(o$inter == FALSE, 1, 2))), ties.method = "min")
  
  # Fit of the final model
  ff <- as.formula(
    paste("Surv(time, status)~",
          ifelse(sum(id == 0) > 0 || sum(sign) > 0,
                 paste(names(data)[attributes(data)$pos$X][which(id %in% c(0, which(sign)))], collapse = "+"),
                 "1")))
  
  fitf <- coxph(ff, data = data)

  res <- coef(fitf)
  names(res) <- names(coef(fitf))
    
  return(list(res, step))
}