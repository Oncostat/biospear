###################################################################################
### Computation of a grid of lambda values from dfmax to 0 penalized coefficients #
###################################################################################

l.Grid <- function(x, y, w, off = 0, alpha = 1, nfold, dfmax, family){
  if(length(off) == 1) off <- rep(0, nrow(x))
  lmax <- glmnet(x = x,
                 y = y,
                 penalty.factor = w,
                 family = family,
                 alpha = alpha,
                 offset = off,
                 standardize = FALSE,
                 dfmax = sum(w == 0))$lambda[1]
  
  lmin <- rev(glmnet(x = x,
                     y = y,
                     penalty.factor = w,
                     family = family,
                     alpha = alpha,
                     offset = off,
                     standardize = FALSE,
                     dfmax = dfmax)$lambda)[1]
  
  lGrid <- c(lmin, lmax)/((1 - (1 / nfold) - 0.1))
  
  return(lGrid)
}