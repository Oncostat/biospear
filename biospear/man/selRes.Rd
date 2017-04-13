\name{selRes}
\alias{selRes}
\title{Evaluation of the selection accuracy of the developed models}
\description{
  This function computes several criteria to assess the selection
accuracy of the developed model. Of note, this function is only
available for simulated datasets for which true biomarkers are known.
}
\usage{
selRes(res)
}
\arguments{
  \item{res}{an object of class '\code{resBMsel}' generated from \code{\link{BMsel}}.}
}
\details{
  Based on the 2x2 contingency table (active vs. inactive / selected vs. unselected),
four selection criteria are provided: \cr
- the false discovery rate (\code{FDR}) that is the proportion of inactive biomarkers
among the selected ones, \cr
- the false non-discovery rate (\code{FNDR}) that is the proportion of active biomarkers
among the unselected ones, \cr
- the false negative rate (\code{FNR}) that is the proportion of unselected
biomarkers among the true active ones,\cr
- and the false positive rate (\code{FPR}) that is the proportion of selected
biomarkers among the inactive ones.\cr
These four criteria are between 0 and 1, and must be minimized. We also provided
two discrimination criteria translating the ability to discard inactive biomarkers
more likely than active ones independently of the tuning parameters: \cr
- the area under the ROC curve (\code{AUC}, depending on the sensitivity and specificity), \cr
- and the area under the precision-recall curve (\code{AUPRC}, depending on the FNR and
FDR). \cr
Of note, the AUPRC is more pertinent than the area under the ROC curve when
there are many more inactive than active biomarkers. These two criteria are between 0 and 1,
and must be maximized.
}
\value{
  \code{matrix} of dimension 6 x the number of implemented methods.
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
  data <- simdata(
    n = 1000, p = 50, q.main = 3, q.inter = 1,
    prob.tt = 0.5, alpha.tt = 0,
    beta.main = c(-0.7, -0.2), beta.inter = -0.5,
    b.corr = 0.6, b.corr.by = 5,
    m0 = 5, wei.shape = 1, recr = 3, fu = 2,
    timefactor = 1)

  res <- BMsel(
    data = data,
    inter = TRUE,
    method = c("alassoL", "stabSel", "lasso-pcvl"),
    folds = 5)

  selAcc <- selRes(res)
}

\keyword{selection}
