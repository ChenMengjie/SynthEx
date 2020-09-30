createSegments <- function(ratio, segmentMethod = c("CBS", "SomaticaEx", "TrendFiltering")){

  if(class(ratio) != "RatioCorrectBiasInTargets" & class(ratio) != "RatioNormalized")
    stop("Invalid input class for createSegments().")

  if(!is.null(ratio$type)){
    if(ratio$type == "WGS" & ratio$Synthetic == TRUE & segmentMethod != "TrendFiltering")
      warning("We recommend to use \'TrendFiltering\' method to segment WGS data not using matched normal.")
  }

  segmentMethod <- match.arg(segmentMethod)
  options(scipen = 50)

  segRes <- .getSegment(ratio$Ratio, 'ratio', segmentMethod)
  if ('normalized' %in% names(ratio$Ratio)) {
    normalizedSeg <- .getSegment(ratio$Ratio, 'normalized', segmentMethod)
  } else {
    normalizedSeg <- NULL
  }

  res <- list(normalizedSeg, segRes, segmentMethod, ratio$Ratio, ratio$Adjust, ratio$TargetBiasStatistics)
  names(res) <- c("segmentNormalized", "segmentUnadjusted", "segmentMethod", "Ratio", "Adjust", "TargetBiasStatistics")
  class(res) <- "Segment"
  return(res)
}
