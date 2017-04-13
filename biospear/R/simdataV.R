#########################################################################
# Simulation function for generating data with biomarkers               #
#########################################################################
# VALIADATION SET                                                       #
#                                                                       #
# - data       the simulated training set                               #
# - NValid     the number of patients of the validation set             #
#########################################################################

simdataV <- function(traindata, Nvalid){  
  
  attr <- attributes(traindata)
  
  dataV <- simdata(
    n = Nvalid,
    p = length(grep("bm", names(traindata))),
    prob.tt = attr$treatment$prob.tt,
    alpha.tt = attr$treatment$alpha.tt,
    q.main = attr$biomarkers$q.main,
    beta.main = attr$biomarkers$beta.main,
    active.main = attr$biomarkers$active.main,
    q.inter = attr$biomarkers$q.inter,
    beta.inter = attr$biomarkers$beta.inter,
    active.inter = attr$biomarkers$active.inter,
    b.corr = attr$biomarkers$correlation,
    b.corr.by = attr$biomarkers$corr.block.size,
    m0 = attr$weibull.parameters$basmed,
    wei.shape = attr$weibull.parameters$shape,
    recr = attr$censoring[1],
    fu = attr$censoring[2],
    timefactor = attr$timefactor
    )
    
  return(dataV)
  
}
#########################################################################
