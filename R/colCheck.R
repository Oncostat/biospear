colCheck <- function(obj, data){
  if(!(is.character(obj) || is.numeric(obj)))
    stop(paste0("The objects 'x', 'y', 'z' or 't' must be defined as character or integer."))
  if(is.character(obj) & length(setdiff(obj, colnames(data))) > 0)
    stop(paste0(paste0(setdiff(obj, colnames(data)), collapse = ", ")," not found in the data set."))
  if(is.numeric(obj)){
    if(sum(obj %% 1) > 0)
      stop(paste0("The objects 'x', 'y', 'z' or 't' must be defined as character or integer."))
    if(min(obj) < 0 || max(obj) > ncol(data))
      stop(paste0(obj[c(which(obj < 0), which(obj > ncol(data)))], ": out of the number of columns in data."))
  }
}
