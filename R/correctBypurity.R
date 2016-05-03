
correctBypurity <- function(segRes, admix){
  segRes$ratio <- 2^segRes[, "log2ratio"]
  corrected <- (segRes$ratio + admix - 1)/admix
  corrected[corrected < 0] <- 0.01
  segRes$correctedRatiobyPurity <- corrected
  segRes$correctedlog2ratiobyPurity <- round(log2(corrected + 0.0001), 4)
  return(segRes)
}


