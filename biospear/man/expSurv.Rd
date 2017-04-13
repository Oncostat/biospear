\name{expSurv}
\alias{expSurv}
\alias{predict.resexpSurv}
\alias{plot.resexpSurv}
\title{Computation of expected survival based on the developed models}
\description{
Based on a prediction model, this function computes expected survival for new patients with associated confidence intervals. The returned object can be plotted to obtain a meaningful graphical visualization.
}
\usage{
expSurv(res, traindata, method, ci.level = .95, boot = FALSE, nboot, smooth,
  pct.group, horizon, trace = TRUE, ncores = 1)

\method{predict}{resexpSurv}(object, newdata, ...)

\method{plot}{resexpSurv}(x, method, pr.group, print.ci, xlim, ylim, xlab, ylab,
...)
}
\arguments{
  \item{res}{an object of class '\code{resBMsel}' generated from \code{\link{BMsel}}}
  \item{traindata}{the \code{data.frame} used for  \code{\link{BMsel}} (training set)}
  \item{method}{method to compute}
  \item{ci.level}{the nominal level for a two-sided confidence interval (CI)}
  \item{boot}{logical value: \code{TRUE} = boostraped CI, \code{FALSE} = analytical CI}
  \item{nboot}{number of bootstrap samples (only used when \code{boot}=\code{TRUE})}
  \item{smooth}{logical value indicating if smoothed B-splines should be computed}
  \item{pct.group}{number or percentile of the prognostic-risk groups. If a single number is provided, all the groups must be defined according to Cox (1957). If percentiles are provided, the sum must be 1 (e.g. 0.164, 0.336, 0.336, 0.164)}
  \item{horizon}{time-horizon to estimate the expected survival probabilities}
  \item{trace}{logical parameter indicating if messages should be printed}
  \item{ncores}{number of PC cores used (useful for the bootstrap CI only)}
  \item{object, x}{an object of class '\code{resexpSurv}' generated from \code{\link{expSurv}}}
  \item{newdata}{\code{\link{data.frame}} containing new patients data. \code{newdata} must have the same variables as \code{traindata}}
  \item{pr.group}{parameter for the \code{plot()} indicating the number of the prognostic-risk group for which the plot will be printed}
  \item{print.ci}{logical parameter for the \code{plot()} indicating if CI will be printed}
  \item{xlim, ylim, xlab, ylab}{usual parameters for \link{plot}}
  \item{...}{other paramaters for \link{predict} or \link{plot}}
}
\details{
   Using an object of class '\code{resBMsel}' from \code{\link{BMsel}}, \code{expSurv} computes expected survival at a given \code{horizon} and constructs confidence intervals thereof either with an analytical (\code{boot} = \code{FALSE}) or non-parametric bootstrap approach (\code{boot} = \code{TRUE}). Smoothed B-splines (logical option \code{smooth}) and categorization of the prognostic score into risk groups (using the option \code{pct.group}) may be used to obtain a meaningful graphical visualization. Predictions for new patients (\code{newdata} data frame) can be computed using \code{predict()}. Graphical visualization can be obtained using \code{plot()}.
}

\value{
  A \code{list} of length three containing the expected survival (\code{surv}) and their corresponding confidence intervals (\code{lower} and \code{upper}).
}
\author{Nils Ternes, Federico Rotolo, and Stefan Michiels\cr
Maintainer: Stefan Michiels \email{stefan.michiels@gustaveroussy.fr}
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

res <- BMsel(
  data = data,
  inter = TRUE,
  method = c("lasso", "lasso-pcvl"))

surv <- expSurv(
  res = res,
  traindata = data,
  boot = TRUE,
  nboot = 20,
  smooth = TRUE,
  pct.group = 4,
  horizon = 4,
  trace = TRUE
)
plot(surv, method = "lasso", pr.group = 3, print.ci = TRUE)

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

  surv <- expSurv(
    res = res,
    traindata = Breast,
    boot = FALSE,
    smooth = TRUE,
    horizon = 4,
    trace = TRUE
  )
  plot(surv, method = "lasso", print.ci = TRUE)

}

########################################
########################################
}


\keyword{survival}
\keyword{estimations}
\keyword{confidence intervals}
