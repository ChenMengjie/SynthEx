createTargetBins <- function(Target, bin.size){
  options(scipen = 50)
  start <- floor(Target[, 2]/bin.size)*bin.size
  end <- ceiling(Target[, 2]/bin.size)*bin.size
  res <- cbind(Target, Target[, 1], start, end)
  all <- NULL
  for(i in 1:nrow(res)){
    if(res[i, 6] - res[i, 5] == bin.size){
      all <- rbind(all, res)
    } else {
      atarget <- seq(res[i, 5], res[i, 6] - bin.size, bin.size)
      len <- length(atarget)
      all <- rbind(all, cbind(res[i, 1:4], atarget, atarget + bin.size))
    }
  }
  all <- as.data.frame(all)
  names(all) <- paste0("V", 1:6)
  return(all)
}












