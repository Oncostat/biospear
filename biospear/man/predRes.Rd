\name{predRes}
\alias{predRes}
\alias{plot.predRes}
\title{Evaluation of the prediction accuracy of the developed models}
\description{
  This function computes several criteria to assess the prediction
  accuracy of the developed model.
}
\usage{
predRes(res, method, traindata, newdata, is.2cv, fold.2cv = 5, horizon,
  trace = TRUE, ncores = 1)

\method{plot}{predRes}(x, method, crit = c("C", "PE", "dC"),
  xlim, ylim, xlab, ylab, col,...)
}
\arguments{
  \item{res}{an object of class '\code{resBMsel}' generated from \code{\link{BMsel}}}
  \item{method}{methods for which prediction criteria are computed}
  \item{traindata}{input \code{\link{data.frame}} used to compute the
\code{res} object}
  \item{newdata}{input \code{\link{data.frame}} not used to compute the
\code{res} object. This object is not mandatory}
  \item{is.2cv}{logical parameter indicating if a double
cross-validation process (2CV) should be performed to mimick an external
validation set}
  \item{fold.2cv}{number of folds for the double cross-validation. Default is 5 and considering a large value should provide extremely large computation time. \code{fold.2cv} must not be considered when \code{is.2cv = FALSE}}
  \item{horizon}{time horizon to compute the prediction criteria}
  \item{trace}{logical parameter indicating if messages should be printed. Default is \code{TRUE}}
  \item{ncores}{number of PC cores used (useful for the double cross-validation only)}
  \item{x}{an object of class '\code{predRes}' generated from \code{\link{predRes}}}
  \item{crit}{parameter indicating the criterion for which the results will be printed (\code{C}: concordance via Uno's C-statistic, \code{PE}: prediction error via integrated Brier score and \code{dC}: delta Uno's C-statistic (for the interaction setting only))}
  \item{xlim, ylim, xlab, ylab, col}{usual parameters for \link{plot}}
  \item{...}{other paramaters for \link{plot}}
}
\details{
To evaluate the accuracy of the selected models, three predictive accuracy measures are implemented: \cr
1) As a measure of overall prediction error of the models, we computed the integrated Brier score. The time-dependent Brier score is a quadratic score based on the predicted time-dependent survival probability (\code{PE}). \cr
2) To evaluate the discrimination of the prediction model, we computed Uno's C-statistic, one of the least biased concordance statistic estimator in the presence of censoring (\code{C}). \cr
3) To evaluate the interaction strength of the models, we computed the absolute difference of the arm-specific Uno's C-statistics (\eqn{\Delta}\code{C}). \cr
For simulated datasets, the predictive accuracy metrics are also computed for the "oracle model" that is the unpenalized Cox proportional hazards model fitted to the truly related biomarkers.
}
\value{
  A \code{list} of the same length of the \code{horizon} considered. Each element of the list contains between 1 and 3 sublists depending on the chosen validation (i.e. training set, internal validation through double cross-validation (2CV) and/or external validation). Each sublist is a \code{matrix} containing the predictive accuracy metrics of the implemented methods.
}
\author{Nils Ternes, Federico Rotolo, and Stefan Michiels\cr
Maintainer: Stefan Michiels \email{stefan.michiels@gustaveroussy.fr}
}
\references{
  Ternes N, Rotolo F and Michiels S.
  Empirical extensions of the lasso penalty to reduce
  the false discovery rate in high-dimensional Cox regression models.
  \emph{Statistics in Medicine} 2016;35(15):2561-73.
  doi:\href{http://onlinelibrary.wiley.com/doi/10.1002/sim.6927/abstract}{10.1002/sim.6927}\cr
  Ternes N, Rotolo F, Heinze G and Michiels S.
  Identification of biomarker-by-treatment interactions in randomized
  clinical trials with survival outcomes and high-dimensional spaces.
  \emph{Biometrical journal}. In press.
  doi:\href{http://onlinelibrary.wiley.com/doi/10.1002/bimj.201500234/abstract}{10.1002/bimj.201500234}
}

\examples{

########################################
# Simulated data set
########################################

data <- simdata(
  n = 1000, p = 50, q.main = 3, q.inter = 2,
  prob.tt = 0.5, alpha.tt = 0,
  beta.main = c(-0.7, -0.2), beta.inter = -0.8,
  b.corr = 0.6, b.corr.by = 5,
  m0 = 5, wei.shape = 1, recr = 3, fu = 2,
  timefactor = 1)

newdata <- simdataV(
  traindata = data,
  Nvalid = 500
)

res <- BMsel(
  data = data,
  inter = TRUE,
  method = c("lasso", "lasso-pcvl"))

summary(res)
summary(res, keep = "xt")

predAcc <- predRes(
  res = res,
  traindata = data,
  newdata = newdata,
  is.2cv = FALSE,
  horizon = 1:4,
  trace = TRUE
)
plot(predAcc, crit = "dC")

########################################
# Real data set
########################################

\dontrun{
  data(Breast)
  dim(Breast)

  pos.tt <- which(names(Breast) %in% "treat")
  pos.y <- rev(which(names(Breast) %in% c("time", "status")))
  pos.x <- (1:ncol(Breast))[-c(pos.tt, pos.y)]

  res <-  BMsel(data = Breast,
                x = pos.x,
                y = pos.y,
                tt = pos.tt,
                inter = FALSE,
                std.x = TRUE,
                folds = 5,
                method = c("lasso", "lasso-pcvl"))

  summary(res)

  predAcc <- predRes(
    res = res,
    traindata = Breast,
    is.2cv = TRUE,
    horizon = 1:4,
    trace = TRUe)
  plot(predAcc, crit = "C")
}

########################################
########################################
}

\keyword{prediction}
