CallWholeGenomeDoubling <- function(segRes, WGD = 1.35, pos.prop.threhold = 0.6, pos.log2ratio.threhold = 0.75){

  segRes$ratio <- 2^segRes[, "log2ratio"]
  subsegRes <- segRes[segRes[, "chr"] %in% c(1:22), ]
  positiveSegments <- subsegRes[subsegRes[, "log2ratio"] > 0.1, ]
  segLens <- (subsegRes[, "end"] - subsegRes[, "start"]+1)
  ploidy <- sum(segLens*subsegRes[, "ratio"])/sum(segLens)
  prop <- sum((positiveSegments[, "end"] - positiveSegments[, "start"])+1)/sum(segLens)

  options(scipen = 50)
    if(ploidy >= WGD | (prop >= pos.prop.threhold & mean(positiveSegments[, "log2ratio"]) >= pos.log2ratio.threhold)){
    WholeGenomeDoubling <- TRUE
  } else {
    WholeGenomeDoubling <- FALSE
  }
  res <- list(WholeGenomeDoubling, ploidy)
  names(res) <- c("WholeGenomeDoubling", "ploidy")
  return(res)

}



