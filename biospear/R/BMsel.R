################################################################################
################################################################################
################################################################################
### Global function for SELecting biomarkers to build a prediction model       #
################################################################################

BMsel <- function(
  ####################################################################
  ######################## *** PARAMETERS *** ########################
  data,             # Dataset
  x,                # Colname(s) or position of the biomarkers
  y,                # Colname(s) or position of the outcome (survival: 1st = time, 2nd = status)
  z,                # Colname(s) or position of the clinical covariates
  tt,               # Colname or position of the treatment
  ### Parametrization
  inter,            # Logical value indicating if interactions should be computed
  std.x,            # Standardization of biomarkers ? (TRUE or FALSE)
  std.i,            # Standardization of interactions ? (TRUE or FALSE)
  std.tt,           # Treatment coding as +/- 0.5 ? (TRUE or FALSE)
  ### Variable selection strategy
  # (*) only available for the predictive model
  method = c(
    'alassoL',      # Adaptive lasso penalty (pre-step weighting via Lasso)
    'alassoR',      # Adaptive lasso penalty (pre-step weighting via Ridge)
    'alassoU',      # Adaptive lasso penalty (pre-step weighting via Univariate)
    'enet',         # Elastic net penalty
    'gboost',       # Gradient boosting
    'glasso',       # Group-lasso penalty (*)
    'lasso',        # Lasso penalty
    'lasso-1se',    # Lasso penalty + One-standard-error rule
    'lasso-AIC',    # Lasso penalty + AIC
    'lasso-BIC',    # Lasso penalty + BIC
    'lasso-HQIC',   # Lasso penalty + HQIC
    'lasso-pct',    # Percentile lasso penalty
    'lasso-pcvl',   # Lasso penalty + pcvl
    'lasso-RIC',    # Lasso penalty + RIC
    'modCov',       # Modified covariates (*)
    'PCAlasso',     # Dimension reduction (PCA) + lasso penalty (*)
    'PLSlasso',     # Dimension reduction (PLS) + lasso penalty (*)
    'ridge',        # Ridge penalty
    'ridgelasso',   # Offset ridge + lasso penalty (*)
    'stabSel',      # Stability selection
    'uniFDR'        # Univariate selection with FDR control
    ),
  folds = 5,        # Number of folds (single value) or a vector of values (size n) assigning observations to folds
  uni.fdr = 0.05,   # Univariate: Threshold FDR to control for multiple testing
  uni.test = 1,     # Univariate: Model comparison approach (1 or 2, 2 only available for predictive setting)
  ss.rando = F,     # Stability selection: Random weights ?
  ss.nsub = 100,    # Stability selection: Number of subsamples
  ss.fsub = 0.5,    # Stability selection: Fraction of samples use in the sampling process
  ss.fwer = 1,      # Stability selection: Parameter to control for the FWER - Number of noise variables
  ss.thr = 0.6,     # Stability selection: threshold used for the stability selection
  dfmax = ncol(data) + 1, # Limit the maximum number of selected variables in the model (useful for large p)
  pct.rep = 1,      # Percentile lasso: Number of replications
  pct.qtl = 0.95,   # Percentile lasso: Percentile of the distribution of lambdas retained
  showWarn = TRUE,  # Show warning messages ?
  trace = TRUE      # Print function's progression
  #####################################################################
  ) {

  #####################################################################
  ### CONTROL CHECKING AND DATA MANIPULATION

  if(is.list(data) == TRUE || is.matrix(data) == TRUE)
    data <- as.data.frame(data)

  arg <- match.call()

  ## Control of the model parametrization
  if(!is.logical(inter))
    stop("The 'inter' parameter must be logical (TRUE or FALSE).")

  ## Control of the selected methods
  method <- match.arg(method, several.ok = TRUE)
  Sel <- method
  int <- intersect(method, c("ridge-lasso", "PCAlasso", "PLSlasso", "glasso", "modCov"))
  if(inter == FALSE && length(int) > 0)
    stop(paste0("\n", paste0(paste0("'",int, "'"), collapse = ", "), " only available for predictive setting (par = 'pred')"))

  ## Control of the data set
  if(missing(data))
    stop("\n The 'data' object is missing.")

  if(is.null(attributes(data)$isSim)){

    # Control of 'x'
    if(missing(x))
      stop("\nThe 'x' object is missing.")
    colCheck(obj = x, data = data)
    if(is.numeric(x))
      x <- colnames(data)[x]
    arg[which(names(arg) == "x")] <- list(x)

    # Control of 'y'
    if(missing(y))
      stop("\nThe 'y' object is missing.")
    colCheck(obj = y, data = data)
    if(is.numeric(y))
      y <- colnames(data)[y]
    arg[which(names(arg) == "y")] <- list(y)

    # Control of 'z', if exists
    if(!missing(z)){
      colCheck(obj = z, data = data)
      if(is.numeric(z))
        z <- colnames(data)[z]
      arg[which(names(arg) == "z")] <- list(z)
    }else{
      z <- NULL
    }

    # Control of 'tt', if exists
    if(inter == TRUE){
      if(missing(tt))
        stop("\nFor a interaction setting ('inter' = TRUE), the treatment object ('tt') must be specified.")
      colCheck(obj = tt, data = data)
      if(is.numeric(tt))
        tt <- colnames(data)[tt]
      arg[which(names(arg) == "tt")] <- list(tt)
    }else{
      if(!missing(tt)){
        colCheck(obj = tt, data = data)
        if(is.numeric(tt))
          tt <- colnames(data)[tt]
        arg[which(names(arg) == "tt")] <- list(tt)
      }else{
        tt <- NULL
      }
    }
    if(sum(duplicated(c(x, y, z, tt))) > 0)
      stop("\nSome duplicates are present in the colnames' objects.")

  }else{

    if(missing(x))
      x <- grep("bm", colnames(data))
    if(is.numeric(x))
      x <- colnames(data)[x]

    y <- c("time", "status")

    if(missing(z)){
      z <- NULL
      if(length(grep("cl", colnames(data))) > 0)
        z <- grep("cl", colnames(data))
    }
    if(is.numeric(z))
      z <- colnames(data)[z]

    if(missing(tt)){
      tt <- NULL
      if(inter == TRUE)
        tt <- "treat"
    }
    if(is.numeric(tt))
      tt <- colnames(data)[tt]
    # if(inter == FALSE){
    #   z <- c(tt, z)
    #   tt <- NULL
    # }
  }

  ## True active biomarkers
  if(!is.null(attributes(data)$isSim))
    effect.tt <- (attributes(data)$treatment$alpha.tt != 0)

  trueBM <- attributes(data)$biomarkers$active.main
  trueBI <- attributes(data)$biomarkers$active.int

  ## Control of the "std.x" parameter
  if(!is.null(attributes(data)$isSim)){
    if(!missing(std.x))
      warning("\nThe 'std.x' parameter will be not taking into account.
              Standardization is already done for generating simulated datasets.")
    std.x <- TRUE
  }else{
    if(missing(std.x)){
      warning("\nThe 'std.x' is missing and is fixed as 'TRUE'.")
      std.x <- TRUE
    }else{
      if(!((std.x) %in% c(TRUE, FALSE)))
        stop("\n The 'std.x' parameter must be either TRUE or FALSE.")
    }
  }

  ## Control of the "std.i" parameter
  if(!missing(std.i)){
    if(!((std.i) %in% c(TRUE, FALSE)))
      stop("\nThe 'std.i' parameter must be either TRUE or FALSE.")
  }
  if(!is.null(attributes(data)$isSim)){
    if(missing(std.i))
      std.i <- FALSE
  }else{
    if(inter == TRUE){
      if(missing(std.i)){
        warning("\nThe 'std.i' is missing and is fixed as 'FALSE'.")
        std.i <- FALSE
      }
    }else{
      if(!missing(std.i))
        warning("\nThe 'std.i' object is not used for the prognostic setting.")
      std.i <- NULL
    }
  }

  ## Control of the "std.tt" parameter
  if(!is.null(attributes(data)$isSim)){
    if(!missing(std.tt))
      warning("\nThe 'std.tt' parameter will be not taking into account.
              Treatment coding is already +/- 0.5 in simulated datasets.")
    std.tt <- TRUE
  }else{
    if(is.null(tt)){
      if(!missing(std.tt))
        warning("\n'std.tt' is not used as 'tt' is missing.")
    }else{
      if(!missing(std.tt)){
        if(!((std.tt) %in% c(TRUE, FALSE)))
          stop("\nThe 'std.tt' parameter must be either TRUE or FALSE.")
      }else{
        if(sum(unique(as.vector(t(data[, tt]))) %in% c(-0.5, +0.5)) != 2){
          warning("\nThe treatment has been recoded as -0.5 and +0.5. To keep the initial coding, please specify std.tt = FALSE")
        }
        std.tt <- TRUE
      }
    }
  }

  ## Control of the 'folds' object
  if(!is.numeric(folds)){
    stop("\nThe 'folds' object must be integer")
  }else{
    if(sum(folds %% 1) > 0)
      stop("\nThe 'folds' object must be integer")
  }
  if(length(folds) == 1){
    if(folds < 3 || folds > nrow(data))
      stop("\nThe number of folds must be between 3 and the number of observations (n).")
    folds <- sample(x = 1:folds, size = nrow(data), replace = TRUE)
  }else{
    if(length(folds) != nrow(data))
      stop("\nThe argument 'folds' must be a unique value or a vector of the size equal to the number of observations (n).")
    folds <- as.numeric(factor(folds))
  }

  list <- list(dfmax = dfmax, uni.fdr = uni.fdr, pct.rep = pct.rep, pct.qtl = pct.qtl, uni.test = uni.test,
               ss.rando = ss.rando, ss.nsub = ss.nsub, ss.fsub = ss.fsub, ss.fwer = ss.fwer, ss.thr = ss.thr)
  i <- lapply(
    X = 1:length(list),
    FUN = function(X){
      if(length(list[[X]]) != 1)
        stop(paste0("\nThe argument '", names(list)[X],"' must be a unique value."))
      })

  if(uni.fdr < 0 || uni.fdr > 1)
    stop("\nThe threshold FDR 'uni.fdr' must be between 0 and 1.")

  if(dfmax < 0)
    stop("\nThe maximum number of selected covariates 'dfmax' must be positive.")

  if(pct.rep < 0)
    stop("\nThe number of replications for the percentile lasso 'pct.rep' must be positive.")

  if(pct.qtl < 0 || pct.qtl > 1)
    stop("\nThe retained percentile of the lambdas 'pct.qtl' must be between 0 and 1.")

  if(uni.test %in% c(1, 2)){
    if(uni.test == 2 && inter == FALSE)
      stop("\nThe argument 'uni.test' cannot be 2 for the prognostic setting (only 1)")
  }else{
    stop("\nThe argument 'uni.test' must be 1 or 2 (2 only available for predictive setting).")
  }

  if(!((ss.rando) %in% c("TRUE", "FALSE")))
    stop("\nThe 'ss.rando' parameter must be either TRUE or FALSE.")

  ss.nsub <- ceiling(ss.nsub)
  if(!(ss.nsub > 10))
    stop("\nThe number of subsamples 'ss.nsub' must be greater than 10.")

  if(ss.fsub < 0 || ss.fsub > 1)
    stop("\nThe fraction of subsamples 'ss.fsub' must be between 0 and 1.")

  if(ss.fwer <= 0)
    stop("\nThe parameter to control the FWER 'ss.fwer' must be positive.")

  if(ss.thr < 0 || ss.thr > 1)
    stop("\nThe threshold 'ss.thr' must be between 0.5 and 1.")

  if(!((showWarn) %in% c(TRUE, FALSE)))
    stop("\nThe 'showWarn' parameter must be either TRUE or FALSE.")
  w <- options()$warn
  if(showWarn == FALSE)
    options(warn = -1)

  if(!((trace) %in% c(TRUE, FALSE)))
    stop("\nThe 'trace' parameter must be either TRUE or FALSE.")

  pos.penExt <- which(method %in% paste0('lasso-', c('1se', 'AIC', 'BIC', 'HQIC', 'pct', 'pcvl', 'RIC')))

  method <- replace(
    x = method,
    list = pos.penExt,
    values = 'penExt')

  o <- list(
    folds = folds,
    uni.fdr = uni.fdr,
    uni.test = uni.test,
    dfmax = dfmax,
    ss.rando = ss.rando,
    ss.nsub = ss.nsub,
    ss.fsub = ss.fsub,
    ss.fwer = ss.fwer,
    ss.thr = ss.thr,
    pct.rep = pct.rep,
    pct.qtl = pct.qtl,
    inter = inter,
    Sel = Sel,
    n.penExt = Sel[pos.penExt],
    isSim = (!is.null(attributes(data)$isSim)))

  #########################################################
  ### Model parametrization

  data <- dataTrans(data = data, tt = tt, x = x, y = y, z = z, std.x = std.x, std.i = std.i, std.tt = std.tt, inter = inter, trace = trace)

  inter <- attributes(data)$inter
  if(!is.null(attributes(data)$na.action))
    o$folds <- o$folds[-c(attributes(data)$na.action)]
  if(o$isSim == FALSE)
    effect.tt <- (inter == TRUE & length(trueBI) > 0)

  ## Oracle model
  internal.trueBM <- internal.trueBI <- NULL
  if(length(trueBM) > 0)
    internal.trueBM <- colnames(data)[attributes(data)$pos$x[which(attributes(data)$inames$x %in% trueBM)]]
  if(length(trueBI) > 0 & inter == TRUE){
    internal.trueBI <- colnames(data)[attributes(data)$pos$xt[which(attributes(data)$inames$x %in% trueBI)]]
  }

    fm <- as.formula(paste0("Surv(time, status) ~ 1",
                          ifelse(effect.tt == TRUE, "+ treat", ""),
                          ifelse(length(internal.trueBM) > 0, paste0(" + ", paste0(internal.trueBM, collapse = " + ")), ""),
                          ifelse(length(internal.trueBI) > 0, paste0(" + ", paste0(internal.trueBI, collapse = " + ")), "")))
  fit.oracle <- coxph(fm, data = data)
  if(!is.null(coef(fit.oracle))){
    oracle <- coef(fit.oracle)
  }

  #########################################################
  ### Biomarker selection

  ### Preliminary step

  if("alassoR" %in% method) method <- c(method, "ridge")
  if(max(c("alassoL", "penExt") %in% method) == 1) method <- c(method, "lasso")

  method <- unique(method)

  ## Matrix specifying sequential order of removing biomarkers
  colStep <- setdiff(c(method, o$n.penExt), "penExt")
  Step <- matrix(NA, nrow = ifelse(o$isSim, sum(attributes(data)$weight != 0), 1), ncol = length(colStep))
  colnames(Step) <- colStep
  if(o$isSim)
    rownames(Step) <- colnames(data)[unlist(attributes(data)$pos[c('x', 'xt')])]

  ### Univariate models for alassoU and/or uniFDR
  if(max(c("alassoU", "uniFDR") %in% method) == 1){
    if(trace == TRUE)
      message(paste0(
        "\rComputing univariate models for: ",
        paste0(c("alassoU", "uniFDR")[c("alassoU", "uniFDR") %in% method], collapse = ", ")), appendLF = TRUE)
    res.unii <- unii(data = data, o = o)
    o$p.uniFDR <- res.unii[[1]]
    o$pr.alsU <- res.unii[[2]]
  }

  ### Computation selection for all methods excepted adaptive lasso and lasso extensions
  pos.sel2 <- which(method %in% c("alassoU", "alassoR", "alassoL", "penExt"))
  if(length(pos.sel2) > 0){
    method2 <- method[pos.sel2]
    method <- method[-pos.sel2]
  }

  res <- res.penExt <- res.als <- list()

  Res <- lapply(method, function(M, o){
    if(trace == TRUE)
      message(paste0(
        "\rComputing selection with method: ", M), appendLF = TRUE)
    flush.console()
    Mres <- eval(parse(text = M))(data = data, o = o)
    return(list(Mres))
  }, o = o)
  names(Res) <- method

  try(o$l.lasso <- Res$lasso[[1]][3], silent = T)
  try(o$pr.alsL <- Res$lasso[[1]][4], silent = T)
  try(o$cv.lasso <- Res$lasso[[1]][5], silent = T)
  try(o$fit.lasso <- Res$lasso[[1]][6], silent = T)
  try(o$pr.alsR <- Res$ridge[[1]][3], silent = T)

  if(length(Res) > 0){
    for(i in 1:length(Res)){
      Step[, method[i]] <- as.numeric(unlist(Res[[i]][[1]][2]))
      res[[i]] <- unlist(Res[[i]][[1]][1])
    }
    names(res) <- method
  }

  ### Computation selection for adaptive lasso and lasso extensions
  if(exists("method2")){

    if("penExt" %in% method2){
      if(trace == TRUE)
        message(paste0(
          "\rComputing selection with method: ", paste0(o$n.penExt, collapse = ", "), " (lasso penalized extensions)"), appendLF = TRUE)
      Res.penExt <- penExt(data = data, o = o)
      Step[, which(colStep %in% o$n.penExt)] <- Step[, 'lasso']
      names(Res.penExt) <- o$n.penExt
      res <- merge.list(res, Res.penExt)
    }

    if(length(grep("alasso", method2)) > 0){

      suf <- gsub("alasso", "", method2[grep("alasso", method2)])

      Res.als <- lapply(suf, function(M, o){

        if(trace == TRUE)
          message(paste0(
            "\rComputing selection with method: alasso", M), appendLF = TRUE)
        flush.console()

        Mres <- alasso(data = data, o = o, pr = M)

        return(Mres)
      }, o = o)

      for(i in 1:length(Res.als)){
        Step[, method2[grep("alasso", method2)][i]] <- as.numeric(unlist(Res.als[[i]][2]))
        res.als[i] <- Res.als[[i]][1]
      }

      names(res.als) <- method2[grep("alasso", method2)]
      res <- merge.list(res, res.als)
    }
  }

  #########################################################
  # Formatting results

  options(warn = w)

  ## Renaming of the biomarkers
  if(exists("oracle"))
    res$oracle <- oracle

  for(i in 1:length(res)){
    names(res[[i]]) <- unlist(
      attributes(data)$inames[c('tt', 'z', 'x', 'xt')], use.names = FALSE)[
        which(colnames(data)[attributes(data)$pos$X] %in% names(res[[i]]))]
  }

  ## Renaming of the Step matrix
  Step <- as.data.frame(Step)
  if(nrow(Step) > 1){
    if(exists("oracle") || o$isSim == TRUE){
      Step[, "oracle"] <- 0
      Step[which(rownames(Step) %in% c(internal.trueBM, internal.trueBI)), "oracle"] <- 1
    }
    rownames(Step) <- unlist(attributes(data)$inames[c('x', 'xt')], use.names = FALSE)
  }

  ## Exclude ridge regression for output
  ridge <- NA
  if(sum(names(res) == "ridge") == 1){
    ridge <- res[[which(names(res) %in% "ridge")]]
    res <- res[-which(names(res) %in% "ridge")]
    Step <- Step[, -which(colnames(Step) %in% "ridge")]
  }

  class(res) <- 'resBMsel'
  attributes(res) <- append(
    x = attributes(res),
    values = list(
      formula = arg,
      isSim = o$isSim,
      inter = inter,
      inames = unlist(attributes(data)$inames[c('tt', 'z', 'x', 'y', 'xt')], use.names = FALSE),
      tnames = attributes(data)$tnames,
      trueBM = trueBM,
      trueBI = trueBI,
      ridge = ridge,
      step = Step,
      std.x = std.x,
      std.i = std.i,
      std.tt = std.tt
      ))
  return(res)
}
################################################################################
################################################################################

print.resBMsel <- function(x, ...) {

  n <- names(x)
  attributes(x) <- NULL
  names(x) <- n
  print(x, ...)

} # end of print

################################################################################

summary.resBMsel <- function(object, show = TRUE, keep = c('tt', 'z', 'x', 'xt'), add.ridge = FALSE, ...) {

  ### Control and checking

  keep <- match.arg(keep, several.ok = TRUE)

  if(!is.logical(show))
    stop("\n'show' must be TRUE or FALSE.")

  if(add.ridge == TRUE){
    if(unique(!is.na(attributes(object)$ridge)) == FALSE){
      warning("\n The 'add.ridge' was not taking into account as the ridge penalty was not computed.")
    }else{
      attr <- attributes(object)
      object <- merge.list(list(ridge = attributes(object)$ridge), object)
      n <- names(object)
      attributes(object) <- attr
      names(object) <- n
    }
  }

  all.names <- unique(unlist(lapply(
    X = 1:length(object),
    FUN = function(X) attributes(object[[X]])$names
  )))
  mat <- matrix(unlist(lapply(
    X = 1:length(object),
    FUN = function(X){
      m <- rep(0, length(all.names))
      names(m) <- all.names
      m[attributes(object[[X]])$names] <- round(object[[X]], 6)
      return(m)
    })), nrow = length(all.names), ncol = length(object), byrow = FALSE)

  rownames(mat) <- all.names
  mat <- merge(
    x = data.frame(
      var = attributes(object)$inames,
      type = as.numeric(factor(attributes(object)$tnames, levels = c("tt", "z", "x", "xt")))),
    y = mat,
    by.x = "var",
    by.y = "row.names",
    all.y = TRUE)
  rownames(mat) <- mat[, "var"]

  mat <- mat[rownames(mat)[order(rank(-rowSums(as.matrix(mat[, -c(1:2), drop = FALSE]) != 0), ties.method = "first"))], ]
  mat <- mat[order(rank(mat[, "type"], ties.method = "first")), ]
  rmat <- rownames(mat)
  mat <- matrix(unlist(c(mat[, -(1:2)])), byrow = FALSE, ncol = length(object), dimnames = list(rmat))
  rownames(mat) <- rmat
  mat <- as.matrix(mat[which(rownames(mat) %in% attributes(object)$inames[which(attributes(object)$tnames %in% keep)]), , drop = FALSE])
  if(nrow(mat) == 0){
    mat <- matrix(rep(0, ncol(mat)), nrow = 1)
    rownames(mat) <- ""
  }

  colnames(mat) <- names(object)

  if(!("oracle" %in% colnames(mat)) & attributes(object)$isSim == TRUE)
    mat <- cbind(mat, oracle = 0)

  if(show == TRUE){
    print(as.table(mat), zero.print = ".")
  }else{
    return(mat)
  }
} # end of summary
################################################################################

