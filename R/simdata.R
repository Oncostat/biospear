#########################################################################
# Simulation function for generating data with biomarkers               #
#########################################################################
#                                                                       #
# - n          the number of patients                                   #
# - p          the number of biomarkers                                 #
# - q.main     the number of active main effects                        #
# - q.int      the number of active interactions                        #
# - t.prob     the treatment assignment probability                     #
# - m0         the baseline median survival time                        #
# - alpha      the treatment effect                                     #
# - beta.main  the effect of the active main effects                    #
# - beta.int   the effect of the active interactions                    #
# - b.corr     the correlation within bm blocks                         #
# - b.corr.by  the size of the blocks of correlated bm's                #
# - wei.shape  the shape parameter of the Weibull distribution          #
# - recr       the recruitment period durations                         #
# - fu         the follow-up period duration                            #
# - timefac    the multiplicative factor for times                      #
#              for instance fix it to 365.25 if m0 is expressed         #
#              in years, but times have to be expressed in days         #
# - act.main   the list of the active main effects (not mandatory)      #
# - act.int    the list of the active interactions (not mandatory)      #
#                                                                       #
#########################################################################

simdata <- function(n, p, q.main, q.inter, prob.tt, m0, alpha.tt, beta.main, beta.inter, b.corr,
                     b.corr.by, wei.shape, recr, fu, timefactor, active.main, active.inter){

  #######################################################################
  ### Parameters, checks and manipulation

  if(missing(q.main)) q.main <- 0
  if(missing(q.inter)) q.inter <- 0
  if(missing(alpha.tt)) alpha.tt <- 0
  if(missing(beta.main)) beta.main <- NULL
  if(missing(beta.inter)) beta.inter <- NULL

  ## Overall control
  list <- list(n = n, p = p, q.main = q.main, q.inter = q.inter, prob.tt = prob.tt, m0 = m0, alpha.tt = alpha.tt,
               beta.main = beta.main, beta.inter = beta.inter, b.corr = b.corr, b.corr.by = b.corr.by,
               wei.shape = wei.shape, recr = recr, fu = fu, timefactor = timefactor)

  invisible(lapply(
    X = 1:length(list),
    FUN = function(X){
      if(!is.null(list[[X]])){
        if(!(is.numeric(list[[X]])))
          stop("\nThe parameter '", names(list)[X],"' must be numerical.")
        if(!(names(list[X]) %in% c("beta.main", "beta.inter"))){
          if(length(list[[X]]) != 1)
            stop("\nThe parameter '", names(list)[X],"' must be a unique value.")
        }
      }
    }))

  ## Sample size and number of biomarkers
  if(n <= 0 || p <= 0)
    stop("\nThe sample size 'n' and the number of biomarkers 'p' must be positive.")

  ## Active biomarkers
  if(q.main < 0 || q.inter < 0)
    stop("\nThe number of active biomarkers 'q.main' and 'q.inter' must be positive or null.")
  if((q.main + q.inter) > p)
    stop("\nThe number of active biomarkers is larger than the overall number of biomarkers 'p'.")
  q.main <- round(q.main, 0)
  q.inter <- round(q.inter, 0)

  ## Effect of active biomarkers
  if(is.null(beta.main)){
    if(q.main > 0)
      stop("\nThe biomarkers' main effects 'beta.main' is missing, with no default.")
  }else{
    if(q.main > 0){
      if(length(beta.main) == 1) beta.main <- rep(beta.main, q.main)
      if(length(beta.main) > 1 && length(beta.main) < q.main){
        if(!length(beta.main) == 2)
          warning(paste0("\nThe length of the parameter 'beta.main' is not 1, 2 or ", q.main, " ('q.main'). ",
                         "Effects are randomly chosen within the range of the given list."))
        beta.main <- runif(n = q.main, min = min(beta.main), max = max(beta.main))
      }
      if(length(beta.main) > q.main)
        stop("\nThe biomarkers main effects 'beta.main' are larger than 'q.main'.")
    }else{
      warning("\nThe biomarkers'main effects 'beta.main' should be not taking into account as 'q.main' = 0.")
    }
  }

  if(is.null(beta.inter)){
    if(q.inter > 0)
      stop("\nThe biomarker-by-treatment interactions 'beta.inter' is missing, with no default.")
  }else{
    if(q.inter > 0){
      if(length(beta.inter) == 1) beta.inter <- rep(beta.inter, q.inter)
      if(length(beta.inter) > 1 && length(beta.inter) < q.inter){
        if(!length(beta.inter) == 2)
          warning(paste0("\nThe length of the parameter 'beta.inter' is not 1, 2 or ", q.inter, " ('q.inter'). ",
                         "Effects are randomly chosen within the range of the given list."))
        beta.inter <- runif(n = q.inter, min = min(beta.inter), max = max(beta.inter))
      }
      if(length(beta.inter) > q.inter)
        stop("\nThe biomarker-by-treatment interactions 'beta.inter' are larger than 'q.inter'.")
    }else{
      warning("\nThe biomarker-by-treatment interactions 'beta.inter' should be not taking into account as 'q.inter' = 0.")
    }
  }

  ## Biomarkers' correlation
  if (!(b.corr >= 0 & b.corr < 1))
    stop("\nThe biomarkers correlation parameter 'b.corr' must be in [0, 1).")
  if (p %% b.corr.by > 0)
    stop(paste0("\nThe size of the blocks of correlated biomarkers ",
                "'b.corr.by' must be a divisor of 'p', ",
                "the total number of biomarkers.\n"))

  ## Treatment assignment probability
  if(prob.tt < 0 || prob.tt > 1)
    stop("\nThe treatment assignment probability 'prob.tt' must be between 0 and 1.")

  ## Baseline median survival time
  if(m0 <= 0)
    stop("\nThe baseline median survival time 'm0' must be positive.")

  ## Weibull shape parameter
  if(wei.shape <= 0)
    stop("\nThe shape parameter 'wei.shape' must be positive.")

  ## Follow-up and recruitment period
  if(fu < 0 || recr < 0)
    stop("\nThe follow-up 'fu' and the recruitment period 'recr' must be positive.")

  ## Multiplicative factor for time
  if(timefactor <= 0)
    stop("\nThe multiplicate factor for time 'timefactor' must be positive.")
  m0 <- m0 * timefactor


  #######################################################################
  ### Data generation

  ## Generation of the treatment assignment
  data <- data.frame(treat = rbinom(n, 1, prob.tt) - .5)

  ## Correlation matrix
  n.blocks <- p %/% b.corr.by
  covMat <- diag(n.blocks) %x% matrix(b.corr^abs(matrix(1:b.corr.by, b.corr.by, b.corr.by, byrow = TRUE) - matrix(1:b.corr.by, b.corr.by, b.corr.by)), b.corr.by, b.corr.by)
  diag(covMat) <- 1

  ## Generation of the biomarkers
  data <- cbind(data, mvrnorm(n, rep(0, p), Sigma = covMat))
  rm(covMat)
  colnames(data)[1 + 1:p] <- paste0("bm", gsub(" ", "0", format(1:p)))

  ## Effects of biomarkers
  HR.main <- HR.inter <- rep(1, nrow(data))
  wact <- wact.main <- wact.inter <- NULL

  if(!missing(active.main)){
    if(q.main != length(active.main))
      stop("The number of biomarkers in 'active.main' is not equal to 'q.main'.")
    if(q.main != length(which(active.main %in% colnames(data))))
      stop("The biomarkers in the 'active.main' object are not found in the data set.")
  }else{
    active.main <- NULL
  }

  if(!missing(active.inter)){
    if(q.inter != length(active.inter))
      stop("The number of biomarkers in 'active.inter' is not equal to 'q.inter'.")
    if(q.inter != length(which(gsub("bi", "bm", active.inter) %in% colnames(data))))
      stop("The biomarkers in 'active.inter' are not all found in the data set.")
  }else{
    active.inter <- NULL
  }

  if(q.main + q.inter != 0){

    if(!is.null(active.main) || !is.null(active.inter)){

      if(q.main > 0 && is.null(active.main)){
        active.inter <- gsub("bi", "bm", active.inter)
        x1 <- (1:p)[-which(colnames(data)[1 + 1:p] %in% active.inter)]
        if(length(x1) > 1)
          x1 <- sample(x = x1, size = q.main, replace = FALSE)
        active.main <- paste0("bm", gsub(" ", "0", format(c(
          p, x1))))[-1]
      }
      if(q.inter > 0 && is.null(active.inter)){
        x2 <- (1:p)[-which(colnames(data)[1 + 1:p] %in% active.main)]
        if(length(x2) > 1)
          x2 <- sample(x = x2, size = q.inter, replace = FALSE)
        active.inter <- paste0("bm", gsub(" ", "0", format(c(
          p, x2))))[-1]
      }
      active.inter <- gsub("bm", "bi", active.inter)
      wact <- c(active.main, gsub("bi", "bm", active.inter))

    }else{

      wact <- paste0("bm", gsub(" ", "0", format(c(
        p, sample(x = 1:p, size = q.main + q.inter, replace = FALSE)))))[-1]
    }

    if(q.main == 0){
      wact.main <- NULL
      wact.inter <- wact
    }else{
      wact.main <- wact[1:q.main]
      HR.main <- apply(
        exp(as.matrix(data[, wact.main, drop=FALSE]) * matrix(beta.main, ncol = q.main, nrow = nrow(data), byrow = T)), 1, prod)
    }

    if(q.inter == 0){
      wact.main <- wact
      wact.inter <- NULL
    }else{
      wact.inter <- wact[q.main + (1 : q.inter)]
      HR.inter[which(data$treat == +0.5)] <- apply(exp(as.matrix(data[which(data$treat == +0.5), wact.inter, drop=FALSE]) * matrix(+beta.inter/2, ncol = q.inter, nrow = nrow(data), byrow = T)[which(data$treat == +0.5),]), 1, prod)
      HR.inter[which(data$treat == -0.5)] <- apply(exp(as.matrix(data[which(data$treat == -0.5), wact.inter, drop=FALSE]) * matrix(-beta.inter/2, ncol = q.inter, nrow = nrow(data), byrow = T)[which(data$treat == -0.5),]), 1, prod)
    }
  }

  # Derived Weibull scale parameters (parameter b = lambda^(-1/rho) in rweibull)
  effect.tt <- alpha.tt
  alpha.tt <- rep(alpha.tt, nrow(data))
  alpha.tt[which(data$treat == +0.5)] <- effect.tt/2
  alpha.tt[which(data$treat == -0.5)] <- -effect.tt/2
  wei.scale <- (m0 * (exp(alpha.tt) * log(2) * HR.main * HR.inter)^ (-1 / wei.shape))

  # Event times
  data$time <- Vectorize(rweibull)(n=1, shape = wei.shape, scale = wei.scale)

  # Censoring times
  cTimes <- fu + recr * runif(n)
  data$status <- as.numeric(data$time <= cTimes)
  data$time <- pmin(data$time, cTimes)

  # Appending attributes
  attributes(data) <- append(attributes(data), list(
    biomarkers = list(
      q.main = q.main,
      active.main = wact.main,
      beta.main = beta.main,
      q.inter = q.inter,
      active.inter = wact.inter,
      beta.inter = beta.inter,
      correlation = b.corr,
      corr.block.size = b.corr.by
      ),
    std = TRUE,
     treatment = list(
      prob.tt = prob.tt,
      alpha.tt = effect.tt
      ),
    weibull.parameters = list(
      shape = wei.shape,
      scale = m0 ^ (-wei.shape) * log(2),
      basmed = m0
      ),
    censoring = c("Recruitment time" = recr,
                  "Minimum follow-up time" = fu),
    timefactor = timefactor,
    isSim = TRUE))

  return(data)
}
#########################################################################
