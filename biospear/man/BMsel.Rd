\name{BMsel}
\alias{BMsel}
\alias{summary.resBMsel}
\title{Biomarker selection in a Cox regression model}
\description{
This function enables to fit a Cox regression model for a prognostic
or a biomarker-by-treatment interaction setting subject to a selection procedure to perform variable selection.
}
\usage{
BMsel(data, x, y, z, tt, inter, std.x, std.i, std.tt, method = c('alassoL',
  'alassoR', 'alassoU', 'enet', 'gboost', 'glasso', 'lasso',
  'lasso-1se', 'lasso-AIC', 'lasso-BIC', 'lasso-HQIC', 'lasso-pct',
  'lasso-pcvl','lasso-RIC', 'modCov', 'PCAlasso', 'PLSlasso', 'ridge',
  'ridgelasso', 'stabSel', 'uniFDR'), folds = 5, uni.fdr = 0.05,
  uni.test = 1, ss.rando = F, ss.nsub = 100, ss.fsub = 0.5, ss.fwer = 1
  , ss.thr = 0.6, dfmax = ncol(data) + 1, pct.rep = 1, pct.qtl = 0.95,
  showWarn = TRUE, trace = TRUE)

\method{summary}{resBMsel}(object, show = TRUE, keep = c('tt', 'z', 'x', 'xt'), add.ridge =
FALSE, ...)
}
\arguments{
  \item{data}{input \code{\link{data.frame}}. Each row is an observation.}
  \item{x}{colnames or position of the biomarkers in \code{data}.}
  \item{y}{colnames or position of the survival outcome in \code{data}. The
first column must be the time and the second must be the indicator.}
  \item{z}{colnames or position of the clinical covariates in \code{data},
  if any.}
  \item{tt}{colnames or position of the treatment in \code{data}, if any.}
  \item{inter}{logical parameter indicating if biomarker-by-treatment
  interactions should be computed.}
  \item{std.x}{logical parameter indicating if the biomarkers should be
standardized (i.e. substracting by the mean and dividing by the standard
deviation of each biomarker).}
  \item{std.i}{logical parameter indicating if the biomarker-by-treatment
interactions should be standardized (i.e. substracting by the mean and
dividing by the standard deviation of each interaction).}
  \item{std.tt}{logical parameter indicating if the treatment should be recoded as +/-0.5.}
  \item{method}{methods computed to perform variable selection and to estimate
the regression coefficients. See the details section to understand all the implemented methods.}
  \item{folds}{number of folds. \code{folds} must be a unique value
comprised between 3 and the sample size (leave-one-out CV, but not
recommended for large datasets) or a vector (same length as the sample size) indicating the fold assignment group of each observation.}
  \item{uni.fdr, uni.test}{specific parameter for the univariate procedure. \code{uni.fdr}: threshold false discovery rate (FDR) to control for multiple testing (Benjamini and Hochberg, 1995), \code{uni.test}: model comparison approach. 1: p-value of the biomarker effect (i.e. main effect for the prognostic setting, or main effect + interaction for the interaction setting), 2: p-value of the interaction (only available for interaction setting).}
  \item{ss.fsub, ss.fwer, ss.nsub, ss.rando, ss.thr}{specific parameter for the stability selection. \code{ss.fsub}: fraction of samples to use in the sampling process (default is 0.5), \code{ss.fwer}: parameter to control for
the family-wise error rate (FWER, i.e. number of noise variables, default is 1), \code{ss.nsub}: number of subsampling (default is 100), \code{ss.rando}: logical parameter indicating if random weights should be added in the lasso penalty (default
is \code{FALSE}), \code{ss.thr}: threshold of the stability probability.
for filtering variable.}
  \item{dfmax}{limit the maximum number of variables in the model. Useful for very large number of covariates to limit the time computation.}
  \item{pct.rep, pct.qtl}{specific parameter for the percentile lasso.
\code{pct.rep}: number of replicates, \code{pct.qtl}: percentile used to estimate the lambda among its generated distribution.}
  \item{showWarn}{logical parameter indicating if warnings should be printed. Default is \code{TRUE}.}
  \item{trace}{logical parameter indicating if messages should be printed. Default is \code{TRUE}.}
  \item{object}{object of class '\code{resBMsel}' return by \code{BMsel}.}
  \item{show}{parameter for the \code{summary()} indicating if the result should be printed. Default is \code{TRUE}.}
  \item{keep}{parameter for the \code{summary()} indicating the type of covariates that should be kept for the summary (\code{tt}: treatment covariate, \code{z}: clinical covariates, \code{x}: biomarker main effects and \code{xt}: biomarker-by-treatment interactions) logical parameter indicating if messages should be printed.}
  \item{add.ridge}{parameter for the \code{summary()} indicating if the ridge penalty should be kept for the summary as no selection is performed. Default is \code{FALSE}.}
  \item{...}{other paramaters for \link{plot} or \link{summary}.}
}
\details{
The \code{method} parameter specifies the chosen implemented approaches. Most of these selection methods are based on the lasso penalty (Tibshirani, 1996). The tuning parameter is usually chosen though the cross-validated log-likelihood criterion (cvl), except for the empirical extensions of the \code{lasso} in which additional penalties to the cvl (given with a suffix, e.g. \code{lasso-pcvl}) were used to estimate the tuning parameter. Other methods based on the lasso were also implemented such as the adaptive lasso (\code{alassoL}, \code{alassoR} and \code{alassoU} for which the last letter indicates the procedure used to estimate the preliminary weights: "\code{L}" for lasso, "\code{R}" for ridge and "\code{U}" for univariate), the elastic-net (\code{enet}) or the stability selection (\code{stabSel}). For the interaction setting, specific methods were implemented: to reduce/control the main effects matrix (i.e. ridge (\code{ridgelasso}) or dimension reduction (\code{PCAlasso} or \code{PLSlasso})), to select or discard main effects and interactions simultaneously (i.e. group-lasso (\code{glasso})), or to include only the interaction part in the model (i.e. \code{modCov}). Some selection methods not based on penalized regression are also proposed: univariate selection (\code{uniFDR}), gradient boosting (\code{gboost}). Finally, even if no selection was performed, the \code{ridge} penalty can be computed. In all cases, clinical covariates were considered as unpenalized or forced in the models.
}
\value{
  An object of class '\code{resBMsel}' containing the list of the selected biomarkers and their estimated regression coefficients for the chosen methods.
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
  beta.main = c(-0.7, -0.2), beta.inter = -0.5,
  b.corr = 0.6, b.corr.by = 5,
  m0 = 5, wei.shape = 1, recr = 3, fu = 2,
  timefactor = 1)

res <- BMsel(
  data = data,
  inter = TRUE,
  method = c("lasso", "lasso-pcvl"))

summary(res)
summary(res, keep = "xt")

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
}

########################################
########################################
}

\keyword{survival}
\keyword{biomarker selection}
\keyword{prognostic}
\keyword{interactions}
