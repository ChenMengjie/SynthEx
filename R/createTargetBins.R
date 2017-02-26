createTargetBins <- function(Target, bin.size){
  options(scipen = 50)
  start <- floor(Target[, 2]/bin.size)*bin.size
  end <- ceiling(Target[, 3]/bin.size)*bin.size
  res <- cbind(Target, Target[, 1], start, end)
  flag <- which(res[, 6] - res[, 5] != bin.size)
  colnames(res) <- paste0("V", 1:6)

  all <- NULL
  for(i in flag){
    atarget <- seq(res[i, 5], res[i, 6] - bin.size, bin.size)
    len <- length(atarget)
    all <- rbind(all, data.frame(res[i, 1:4], atarget, atarget + bin.size, row.names = NULL))
  }
  colnames(all) <- paste0("V", 1:6)

  final <- rbind(res[-flag, ], all)
  return(final)
}












