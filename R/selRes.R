################################################################################
################################################################################
################################################################################
### Selection accuracy of the methods                                          #
### WARNING : Only available for simulation datasets (known gold-standard)     #
################################################################################

selRes <- function(
  ####################################################################
  ######################## *** PARAMETERS *** ########################
  res              # Object of class 'resBMsel'
  ######################## *** PARAMETERS *** ########################  
){

  ####################################################################
  ### CHECKING AND MANIPULATION
  
  if(attributes(res)$isSim == FALSE)
    stop("\nOnly available for the results of a simulated data set.")

  ## 'Res' matrix for the computation of FDR and FNR
  Res <- data.frame(summary(res, show = FALSE, keep = ifelse(attributes(res)$inter == FALSE, "x", "xt")))

  ## 'Step' matrix for the computation of AUC and AUPRC
  Step <- data.frame(attributes(res)$step)
  if(attributes(res)$inter == TRUE)
    Step <- Step[-(1:(nrow(Step)/2)), ]
  Step <- Step[, colnames(Res)]

  nmeth <- ncol(Res) - 1
  allbiom <- attributes(res)$inames[which(attributes(res)$tnames == ifelse(attributes(res)$inter == FALSE, "x", "xt"))]

  ####################################################################
  ### COMPUTATION OF SELECTION CRITERIA
  
  ## FALSE DISCOVERY RATE
  FDR <- unlist(lapply(
    X = 1:nmeth,
    function(X) length(setdiff(rownames(Res)[which(Res[, X] != 0)], rownames(Res)[which(Res[, 'oracle'] != 0)])) / length(which(Res[, X] != 0))
  ))
  FDR <- ifelse(is.na(FDR), 0, FDR)

  ## FALSE NON-DISCOVERY RATE
  FNDR <- unlist(lapply(
    X = 1:nmeth,
    function(X) 1 - (length(setdiff(setdiff(allbiom, rownames(Res)[which(Res[, X] != 0)]), rownames(Res)[which(Res[, 'oracle'] != 0)])) / (length(allbiom) - length(which(Res[, X] != 0))))
  ))
  FNDR <- ifelse(is.na(FNDR), 0, FNDR)

  ## FALSE NEGATIVE RATE
  FNR <- rep(NA, nmeth)
  if(sum(Res[, "oracle"]) != 0){
    FNR <- unlist(lapply(
      X = 1:nmeth,
      function(X) 1 - (length(intersect(rownames(Res)[which(Res[, 'oracle'] != 0)], rownames(Res)[which(Res[, X] != 0)])) / sum(Res[,"oracle"] != 0))
    ))
  }

  ## FALSE POSITIVE RATE
  FPR <- rep(NA, nmeth)
  if(sum(Res[, "oracle"]) != nrow(Res)){
    FPR <- unlist(lapply(
      X = 1:nmeth,
      function(X) 1 - (length(setdiff(setdiff(allbiom, rownames(Res)[which(Res[, 'oracle'] != 0)]), rownames(Res)[which(Res[, X] != 0)])) / (length(allbiom) - sum(Res[,"oracle"] != 0)))
    ))
  }

  ## AREA UNDER THE ROC CURVE
  AUC <- rep(NA, nmeth)
  if(sum(Step[, "oracle"]) != 0){
    AUC <- unlist(lapply(
      X = 1:nmeth,
      function(X) auc(roc(response = Step[, 'oracle'], predictor = Step[, X]))
    ))
  }

  ## AREA UNDER THE PRECISION-RECALL CURVE (Davis and Goadrich, 2006)
  AUPRC <- rep(NA, nmeth)
  if(sum(Step[, "oracle"]) != 0){
    AUPRC <- unlist(lapply(
      X = 1:nmeth,
      function(X) pr.curve(scores.class0 = Step[, X], weights.class0 = Step[, 'oracle'])$auc.integral
    ))
  }

  selRes <- round(rbind(FDR, FNDR, FNR, FPR, AUC, AUPRC), 4)
  colnames(selRes) <- names(res)[1:nmeth]
  return(selRes)
}
