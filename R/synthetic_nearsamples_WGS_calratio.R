synthetic_nearsamples_WGS_calratio <- function(tumor, counts, bin.size = 100000, rm.centromere = TRUE,
                                   centromereBins = NULL, K = 5, reads.threshold = 50, chrX = FALSE){

  options(scipen = 50)
  sampleData <- tumor
  if(substr(sampleData[1, 1], 1, 3) == "chr") {
    sampleData[, 1] <- gsub("chr", "", sampleData[, 1])
  }
  sampleData[, 1] <- gsub("X", "23", sampleData[, 1])
  if(nrow(counts) != nrow(sampleData)){
    stop("Input data and \"counts\" size don't match!")
  }

  all.value <- NULL
  for(i in 1:ncol(counts)){
    normal <- counts[, i]
    sampleData[sampleData[, "reads"] < reads.threshold, "reads"] <- 0
    normal[normal < reads.threshold] <- 0
    ratio <- sampleData[, "reads"]/normal
    ratio <- ratio[is.finite(ratio) & ratio != 0 ]
    ratio <- ratio/median(ratio, na.rm = T)
    log.ratio <- log2(ratio[!is.na(ratio)]+0.001)
    diff.sumsquare <- abs(diff(log.ratio))^2
    variance <- mean(diff.sumsquare[diff.sumsquare < quantile(diff.sumsquare, 0.9)])
    all.value <- c(all.value, variance)
  }
  minIs <- order(all.value)[1:K]
  normal <- apply(counts[, minIs], 2, median)
  normal[normal < reads.threshold] <- 0
  ratio <- sampleData[, "reads"]/normal
  ratio <- ratio/median(ratio[is.finite(ratio) & ratio != 0], na.rm = T)
  ratio[is.infinite(ratio) | is.nan(ratio)] <- NA
  ratio.res <- data.frame(sampleData[, c(1:3)], ratio)

  if(rm.centromere == TRUE) {
    if(is.null(centromereBins)){
      if(!bin.size %in% c(10000, 25000, 50000, 100000)){
        stop(paste0("SynthEx doesn't have centromere bins pre-calculated for bin size of", bin.size, "; Please use createCentromereBins()
                    to generate the required file or consider to use another bin size."))
      } else {
        data(CentromereAnnotations)
        ss <- paste0("centromere <- CentromereAnnotations$bin", bin.size)
        eval(parse(text=ss))
      }
    } else {
      centromere <- read.delim(centromereBins, header = F, as.is = T)
    }
    centromere.IDs <- paste0(centromere[, 1], ":", centromere[, 2])
    ratio.IDs <- paste0(ratio.res[, "chr"], ":", ratio.res[, "start"])
    ratio.res <- ratio.res[! ratio.IDs%in%centromere.IDs, ]
  }

  if(chrX == FALSE){
    ratio.res <- ratio.res[ratio.res[, "chr"] != 23 & ratio.res[, "chr"] != 24, ]
  }

  res <- list(ratio.res, TRUE)
  names(res) <- c("Ratio", "WGS")
  class(res) <- "WGSRatio"
  return(res)

}
