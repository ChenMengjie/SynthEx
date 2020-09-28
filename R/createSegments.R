.getSegment <- function(df, value_field = c("ratio", "normalized"), segmentMethod = c("CBS", "SomaticaEx", "TrendFiltering")) {

  subRatio <- df[!is.na(df[, value_field]) & is.finite(df[, value_field]) &
                      df[, value_field] > 0.1, c("chr", "start", "end", value_field)]
  bin.size <- df$end[1] - df$start[1] + 1
  if(segmentMethod == "CBS"){
    suppressPackageStartupMessages(require(DNAcopy))
    seg <- segment(CNA(log2(subRatio[, value_field]+0.0001), subRatio[, "chr"], subRatio[, "start"]), verbose = 0)
    segRes <- seg$output[, c(2:4, 6)]
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] + bin.size - 1
  }

  if(segmentMethod == "SomaticaEx"){
    segRes <- larsCBSsegment(cbind(subRatio[, c("chr", "start")], log2(subRatio[, value_field]+0.0001)), verbose = F, k = 10)
    segRes[, 4] <- round(as.numeric(as.character(segRes[, 4])), 4)
    segRes[, 2:3] <- apply(segRes[, 2:3], 2, as.numeric)
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] - 1
    segRes[, 1] <- as.character(segRes[, 1])
  }

  if(segmentMethod == "TrendFiltering"){
    suppressPackageStartupMessages(require(flsa))
    lambda2 <- 0:60/10
    ratioDiff <- diff(subRatio[, value_field])
    ratioDiff <- ratioDiff[!is.nan(ratioDiff) & is.finite(ratioDiff)]
    noiseLevel <- sqrt(sum(abs(ratioDiff)^2)/length(ratioDiff))
    flsaRes <- flsa(log2(subRatio[, value_field]+0.0001), lambda2 = lambda2)
    sdList <- apply(flsaRes, 1, function(x){
      return(sd(x - log2(subRatio[, value_field]+0.0001)))
    })
    best.num <- which.min(abs(sdList - noiseLevel))
    denoised.ratio <- flsaRes[best.num, ]
    segRes <- .chr.tosegment(round(denoised.ratio, 2), subRatio)
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] - 1
  }
  segRes <- data.frame(segRes)

  colnames(segRes) <-  c("chr", "start", "end", "log2ratio")
  segRes <- segRes[segRes[, "end"] != 0, ]
  segRes <- segRes[segRes[, "end"] > segRes[, "start"], ]
  return(segRes)
}

createSegments <- function(ratio, segmentMethod = c("CBS", "SomaticaEx", "TrendFiltering")){

  if(class(ratio) != "RatioCorrectBiasInTargets" & class(ratio) != "RatioNormalized")
    stop("Invalid input class for createSegments().")

  if(!is.null(ratio$type)){
    if(ratio$type == "WGS" & ratio$Synthetic == TRUE & segmentMethod != "TrendFiltering")
      warning("We recommend to use \'TrendFiltering\' method to segment WGS data not using matched normal.")
  }

  segmentMethod <- match.arg(segmentMethod)
  options(scipen = 50)

  ratioNormalized <- ratio$Ratio
  segRes <- .getSegment(ratio$Ratio, 'ratio', segmentMethod)
  if ('normalized' %in% ratio$Ratio) {
    normalizedSeg <- .getSegment(ratio$Ratio, 'normalized', segmentMethod)
  } else {
    normalizedSeg <- NULL
  }

  res <- list(normalizedSeg, segRes, segmentMethod, ratio$Ratio, ratio$Adjust, ratio$TargetBiasStatistics)
  names(res) <- c("segmentNormalized", "segmentUnadjusted", "segmentMethod", "Ratio", "Adjust", "TargetBiasStatistics")
  class(res) <- "Segment"
  return(res)
}


