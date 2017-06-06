################################################################################
### Internal function for univariate models                                    #
################################################################################

unii <- function(data, o){
  
  if(o$inter == FALSE) id <- c(rep(0, sum(attributes(data)$weights == 0)), 1:sum(attributes(data)$weights == 1))    
  if(o$inter == TRUE) id <- c(rep(0, sum(attributes(data)$weights == 0)), rep(1:(sum(attributes(data)$weights == 1) / 2), 2))   
  
  f0 <- as.formula(
    paste("Surv(time, status)~",
          ifelse(sum(id == 0) > 0,
                 paste(names(data)[attributes(data)$pos$X][which(id == 0)], collapse = "+"),
                 "1"))
  )
  
  # Computation of (p) p-values
  res <- unlist(lapply( 
    X = 1:max(id),
    FUN = function(X){       
      
      if(o$uni.test == 2)
        f0 <- as.formula(paste("Surv(time, status)~", 
                               if(min(id) == 0) paste(paste(names(data)[attributes(data)$pos$X][which(id == 0)], collapse = "+"), "+") ,
                               paste(names(data)[attributes(data)$pos$X][which(id %in% X)[1]], collapse = "+")))
      
      f <- as.formula(paste("Surv(time, status)~", 
                            if(min(id) == 0) paste(paste(names(data)[attributes(data)$pos$X][which(id == 0)], collapse = "+"), "+") ,
                            paste(names(data)[attributes(data)$pos$X][which(id %in% X)], collapse = "+")))

      fit0 <- coxph(f0, data = data)
      fit <- coxph(f, data = data)
      p <- anova(fit0, fit, test = "F")$P[2]
      coef <- coef(fit)[c(grep("bm", names(coef(fit))), grep("bi", names(coef(fit))))]

      return(c(coef, p))  
    }))
  p.uniFDR <- res[(names(res) == "")]
  pr.alsU <- res[!(names(res) == "")]
  pr.alsU <- pr.alsU[c(grep("bm", names(pr.alsU)), grep("bi", names(pr.alsU)))]
  
  return(list(p.uniFDR, pr.alsU))
}