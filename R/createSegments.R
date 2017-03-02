createSegments <- function(ratio, segmentMethod = c("CBS", "SomaticaEx", "TrendFiltering")){

  if(class(ratio) != "RatioNormalized")
    stop("Invalid input class for createSegments().")

  if(!is.null(ratio$type)){
    if(ratio$type == "WGS" & ratio$Synthetic == TRUE & segmentMethod != "TrendFiltering")
      warning("We recommend to use \'TrendFiltering\' method to segment WGS data not using matched normal.")
  }

  segmentMethod <- match.arg(segmentMethod)
  options(scipen = 50)

  ratioNormalized <- ratio$Ratio
  subratioNormalized <- ratioNormalized[!is.na(ratioNormalized[, "normalized"]) & is.finite(ratioNormalized[, "normalized"]) &
                                             ratioNormalized[, "normalized"] > 0.05, c("chr", "start", "end", "normalized")]

  bin.size <- ratioNormalized$end[1] - ratioNormalized$start[1] + 1
  if(segmentMethod == "CBS"){
    suppressPackageStartupMessages(require(DNAcopy))
    #seg <- segment(CNA(subratioNormalized[, "normalized"], subratioNormalized[, "chr"], subratioNormalized[, "start"]), verbose = 0)
    #segRes <- seg$output[, c(2:4, 6)]
    #segRes[, 3] <- segRes[, 3] + bin.size - 1
    #segRes[, 4] <- round(log2(segRes[, 4]+0.0001), 4)
    seg <- segment(CNA(log2(subratioNormalized[, "normalized"]+0.0001), subratioNormalized[, "chr"], subratioNormalized[, "start"]), verbose = 0)
    segRes <- seg$output[, c(2:4, 6)]
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] + bin.size - 1
  }

  if(segmentMethod == "SomaticaEx"){

    #segRes <- larsCBSsegment(subratioNormalized[, c("chr", "start", "normalized")], verbose = F, k = smoothk)
    #segRes[, 4] <- round(log2(as.numeric(as.character(segRes[, 4]))+0.0001), 4)

    segRes <- larsCBSsegment(cbind(subratioNormalized[, c("chr", "start")], log2(subratioNormalized[, "normalized"]+0.0001)), verbose = F, k = 10)
    segRes[, 4] <- round(as.numeric(as.character(segRes[, 4])), 4)
    segRes[, 2:3] <- apply(segRes[, 2:3], 2, as.numeric)
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] - 1
    segRes[, 1] <- as.character(segRes[, 1])
  }

  if(segmentMethod == "TrendFiltering"){
    suppressPackageStartupMessages(require(flsa))
    lambda2 <- 0:60/10
    ratioDiff <- diff(subratioNormalized[, "normalized"])
    ratioDiff <- ratioDiff[!is.nan(ratioDiff) & is.finite(ratioDiff)]
    noiseLevel <- sqrt(sum(abs(ratioDiff)^2)/length(ratioDiff))
    #flsaRes <- flsa(subratioNormalized[, "normalized"], lambda2 = lambda2)
    #sdList <- apply(flsaRes, 1, function(x){
    #  return(sd(x - subratioNormalized[, "normalized"]))
    #})
    flsaRes <- flsa(log2(subratioNormalized[, "normalized"]+0.0001), lambda2 = lambda2)
    sdList <- apply(flsaRes, 1, function(x){
      return(sd(x - log2(subratioNormalized[, "normalized"]+0.0001)))
    })
    best.num <- which.min(abs(sdList - noiseLevel))
    denoised.ratio <- flsaRes[best.num, ]
    segRes <- .chr.tosegment(round(denoised.ratio, 2), subratioNormalized)
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] - 1
    #segRes[, 4] <- round(log2(segRes[, 4]+0.0001), 4)
  }

  normalizedSeg <- data.frame(segRes)

  subRatio <- ratioNormalized[!is.na(ratioNormalized[, "ratio"]) & is.finite(ratioNormalized[, "ratio"]) &
                                  ratioNormalized[, "ratio"] > 0.1, c("chr", "start", "end", "ratio")]

  if(segmentMethod == "CBS"){
   # seg <- segment(CNA(subRatio[, "ratio"], subRatio[, "chr"], subRatio[, "start"]), verbose = 0)
   # segRes <- seg$output[, c(2:4, 6)]
   # segRes[, 3] <- segRes[, 3] + bin.size - 1
   # segRes[, 4] <- round(log2(segRes[, 4]+0.0001), 4)
    seg <- segment(CNA(log2(subRatio[, "ratio"]+0.0001), subRatio[, "chr"], subRatio[, "start"]), verbose = 0)
    segRes <- seg$output[, c(2:4, 6)]
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] + bin.size - 1
  }

  if(segmentMethod == "SomaticaEx"){
   # segRes <- larsCBSsegment(subRatio[, c("chr", "start", "ratio")], verbose = F, k = 10)
   # segRes[, 4] <- round(log2(as.numeric(as.character(segRes[, 4]))+0.0001), 4)
    segRes <- larsCBSsegment(cbind(subRatio[, c("chr", "start")], log2(subRatio[, "ratio"]+0.0001)), verbose = F, k = 10)
    segRes[, 4] <- round(as.numeric(as.character(segRes[, 4])), 4)
    segRes[, 2:3] <- apply(segRes[, 2:3], 2, as.numeric)
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] - 1
    segRes[, 1] <- as.character(segRes[, 1])
  }

  if(segmentMethod == "TrendFiltering"){
    suppressPackageStartupMessages(require(flsa))
    lambda2 <- 0:60/10
    ratioDiff <- diff(subRatio[, "ratio"])
    ratioDiff <- ratioDiff[!is.nan(ratioDiff) & is.finite(ratioDiff)]
    noiseLevel <- sqrt(sum(abs(ratioDiff)^2)/length(ratioDiff))
    #flsaRes <- flsa(subRatio[, "ratio"], lambda2 = lambda2)
    #sdList <- apply(flsaRes, 1, function(x){
    #  return(sd(x - subRatio[, "ratio"]))
    #})
    flsaRes <- flsa(log2(subRatio[, "ratio"]+0.0001), lambda2 = lambda2)
    sdList <- apply(flsaRes, 1, function(x){
      return(sd(x - log2(subRatio[, "ratio"]+0.0001)))
    })
    best.num <- which.min(abs(sdList - noiseLevel))
    denoised.ratio <- flsaRes[best.num, ]
    segRes <- .chr.tosegment(round(denoised.ratio, 2), subRatio)
    segRes[, 2] <- segRes[, 2] - 1
    segRes[, 3] <- segRes[, 3] - 1
    segRes[, 4] <- round(log2(segRes[, 4]+0.0001), 4)
  }
  segRes <- data.frame(segRes)

  colnames(segRes) <- colnames(normalizedSeg) <- c("chr", "start", "end", "log2ratio")

  segRes <- segRes[segRes[, "end"] != 0, ]
  normalizedSeg <- normalizedSeg[normalizedSeg[, "end"] !=0, ]
  segRes <- segRes[segRes[, "end"] > segRes[, "start"], ]
  normalizedSeg <- normalizedSeg[normalizedSeg[, "end"]  > normalizedSeg[, "start"], ]

  res <- list(normalizedSeg, segRes, segmentMethod, ratio$Ratio, ratio$Adjust, ratio$TargetBiasStatistics)
  names(res) <- c("segmentNormalized", "segmentUnadjusted", "segmentMethod", "Ratio", "Adjust", "TargetBiasStatistics")
  class(res) <- "Segment"
  return(res)
}


