[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/biospear)](https://cran.r-project.org/package=biospear)

# biospear
Biomarker selection in penalized regression models

## Aim
We proposed biospear, a useful R tool for developing and validating prediction models, estimate expected survival of patients and visualize them graphically. A function generating simulated survival data set is also provided.

# Install directly from github
#Could be needed to install some of these packages first
##install.packages(c("cobs", "grplasso", "mboost", "plsRcox", "pROC", "PRROC", "RCurl", "survAUC"))

#install.packages(c("dplyr", "glmnet", "corpcor", "survival","statnet.common","devtools"))

source("http://bioconductor.org/biocLite.R")

biocLite("survcomp")

library(devtools)

install_github("oncostat/biospear",type="source")


## References
Ternès N, Rotolo F, Michiels S. biospear: an R package for biomarker selection in penalized Cox regression. Bioinformatics. 2018 Jan 1;34(1):112-113. doi: [10.1093/bioinformatics/btx560] (http://dx.doi.org/10.1093/bioinformatics/btx560).

Ternes N, Rotolo F and Michiels S.
Empirical extensions of the lasso penalty to reduce 
the false discovery rate in high-dimensional Cox regression models.
*\emph{*Statistics in Medicine* 2016;**35**(15):2561-2573.
DOI: [10.1002/sim.6927](http://dx.doi.org/10.1002/sim.6927)

Ternes N, Rotolo F, Heinze G and Michiels S.
Identification of biomarker-by-treatment interactions in randomized
clinical trials with survival outcomes and high-dimensional spaces.
*Biometrical journal* 2017 **59**(4):685-701.
DOI: [10.1002/bimj.201500234](http://dx.doi.org/10.1002/bimj.201500234)

Ternes N, Rotolo F and Michiels S.
Robust estimation of the expected survival probabilities from high-dimensional Cox models with biomarker-by-treatment interactions in randomized clinical trials.
*BMC Medical Research Methodology*. 2017 **17**:83.
DOI: [10.1002/bimj.201500234](http://dx.doi.org/10.1186/s12874-017-0354-0)
