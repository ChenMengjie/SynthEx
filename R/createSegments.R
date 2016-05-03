createSegments <- function(ratio, segmentMethod = c("CBS", "SomaticaEx", "TrendFiltering"),
                           smoothk = 10, ratio.cutoff = 0.05){

  if(class(ratio) != "RatioNormalized" & class(ratio) != "RatioCorrectBiasInTargets" & class(ratio) != "WGSRatio")
    stop("Invalid input class for createSegments().
         The input of createSegments() should be the output of correctBias() or normalization().")

  if(!is.null(ratio$WGS)){
    if(ratio$WGS == TRUE & segmentMethod != "TrendFiltering")
      warning("We recommend to test \'TrendFiltering\' method to on WGS data.")
  }

  segmentMethod <- match.arg(segmentMethod)
  options(scipen = 50)

  ratioNormalized <- ratio$Ratio
  bin.size <- ratioNormalized$end[1] - ratioNormalized$start[1] + 1

  if(class(ratio) == "RatioNormalized" ){
    subratioNormalized <- ratioNormalized[!is.na(ratioNormalized[, "normalized"]) & is.finite(ratioNormalized[, "normalized"]) &
                                             ratioNormalized[, "normalized"] >= ratio.cutoff, c("chr", "start", "end", "normalized")]

    if(segmentMethod == "CBS"){
      require(DNAcopy)
      seg <- DNAcopy::segment(CNA(subratioNormalized[, "normalized"], subratioNormalized[, "chr"], subratioNormalized[, "start"]), verbose = 0)
      segRes <- seg$output[, c(2:4, 6)]
      segRes[, 3] <- segRes[, 3] + bin.size - 1
      segRes[, 4] <- round(log2(segRes[, 4]+0.0001), 4)
    }

    if(segmentMethod == "SomaticaEx"){
      segRes <- larsCBSsegment(subratioNormalized[, c("chr", "start", "normalized")], verbose = F, k = smoothk)
      segRes[, 4] <- round(log2(as.numeric(as.character(segRes[, 4]))+0.0001), 4)
      segRes[, 2:3] <- apply(segRes[, 2:3], 2, as.numeric)
      segRes[, 3] <- segRes[, 3] - 1
      segRes[, 1] <- as.character(segRes[, 1])
    }

    if(segmentMethod == "TrendFiltering"){
      require(flsa)
      lambda2 <- 0:60/10
      ratioDiff <- diff(subratioNormalized[, "normalized"])
      ratioDiff <- ratioDiff[!is.nan(ratioDiff) & is.finite(ratioDiff)]
      noiseLevel <- sqrt(sum(abs(ratioDiff)^2)/length(ratioDiff))
      flsaRes <- flsa(subratioNormalized[, "normalized"], lambda2 = lambda2)
      sdList <- apply(flsaRes, 1, function(x){
        return(sd(x - subratioNormalized[, "normalized"]))
      })
      best.num <- which.min(abs(sdList - noiseLevel))
      denoised.ratio <- flsaRes[best.num, ]
      segRes <- .chr.tosegment(round(denoised.ratio, 2), subratioNormalized)
      segRes[, 3] <- segRes[, 3] - 1
      segRes[, 4] <- round(log2(segRes[, 4]+0.0001), 4)
    }
    normalizedSeg <- data.frame(segRes)
    colnames(normalizedSeg) <- c("chr", "start", "end", "log2ratio")
  }

  subRatio <- ratioNormalized[!is.na(ratioNormalized[, "ratio"]) & is.finite(ratioNormalized[, "ratio"]) &
                                  ratioNormalized[, "ratio"] >= ratio.cutoff, c("chr", "start", "end", "ratio")]

  if(segmentMethod == "CBS"){
    seg <- segment(CNA(subRatio[, "ratio"], subRatio[, "chr"], subRatio[, "start"]), verbose = 0)
    segRes <- seg$output[, c(2:4, 6)]
    segRes[, 3] <- segRes[, 3] + bin.size - 1
    segRes[, 4] <- round(log2(segRes[, 4]+0.0001), 4)
  }

  if(segmentMethod == "SomaticaEx"){
    segRes <- larsCBSsegment(subRatio[, c("chr", "start", "ratio")], verbose = F, k = 10)
    segRes[, 4] <- round(log2(as.numeric(as.character(segRes[, 4]))+0.0001), 4)
    segRes[, 2:3] <- apply(segRes[, 2:3], 2, as.numeric)
    segRes[, 3] <- segRes[, 3] - 1
    segRes[, 1] <- as.character(segRes[, 1])
  }

  if(segmentMethod == "TrendFiltering"){
    require(flsa)
    lambda2 <- 0:60/10
    ratioDiff <- diff(subRatio[, "ratio"])
    ratioDiff <- ratioDiff[!is.nan(ratioDiff) & is.finite(ratioDiff)]
    noiseLevel <- sqrt(sum(abs(ratioDiff)^2)/length(ratioDiff))
    flsaRes <- flsa(subRatio[, "ratio"], lambda2 = lambda2)
    sdList <- apply(flsaRes, 1, function(x){
      return(sd(x - subRatio[, "ratio"]))
    })
    best.num <- which.min(abs(sdList - noiseLevel))
    denoised.ratio <- flsaRes[best.num, ]
    segRes <- .chr.tosegment(round(denoised.ratio, 2), subRatio)
    segRes[, 3] <- segRes[, 3] - 1
    segRes[, 4] <- round(log2(segRes[, 4]+0.0001), 4)
  }
  segRes <- data.frame(segRes)

  colnames(segRes) <- c("chr", "start", "end", "log2ratio")

  #segRes <- segRes[segRes[, "end"] != 0, ]
  #normalizedSeg <- normalizedSeg[normalizedSeg[, "end"] !=0, ]
  #segRes <- segRes[segRes[, "end"] > segRes[, "start"], ]
  #normalizedSeg <- normalizedSeg[normalizedSeg[, "end"]  > normalizedSeg[, "start"], ]

  if(class(ratio) == "RatioNormalized" ){
    if(!is.null(ratio$WGS)){
      res <- list(normalizedSeg, segRes, segmentMethod, ratio$Ratio, ratio$Adjust, ratio$WGS)
      names(res) <- c("segmentNormalized", "segmentUnadjusted", "segmentMethod", "Ratio", "Adjust", "WGS")
    } else {
      res <- list(normalizedSeg, segRes, segmentMethod, ratio$Ratio, ratio$Adjust, ratio$TargetBiasStatistics)
      names(res) <- c("segmentNormalized", "segmentUnadjusted", "segmentMethod", "Ratio", "Adjust", "TargetBiasStatistics")
    }
  } else if (class(ratio) == "RatioCorrectBiasInTargets") {
    res <- list(segRes, segmentMethod, ratio$Ratio, ratio$TargetBiasStatistics)
    names(res) <- c("segmentUnadjusted", "segmentMethod", "Ratio", "TargetBiasStatistics")
  } else {
    res <- list(segRes, segmentMethod, ratio$Ratio, ratio$WGS)
    names(res) <- c("segmentUnadjusted", "segmentMethod", "Ratio", "WGS")
  }

  class(res) <- "Segment"
  return(res)
}



