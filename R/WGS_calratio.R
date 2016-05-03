WGS_calratio <- function(tumor, normal, bin.size = 100000, rm.centromere = TRUE,
                 centromereBins = NULL, reads.threshold = 50, chrX = FALSE){

  options(scipen = 50)

  x <- tumor
  y <- normal
  x[x[, "reads"] <= reads.threshold, "reads"] <- 0
  y[y[, "reads"] <= reads.threshold, "reads"] <- 0

  factor <- sum(x[, "reads"])/sum(y[, "reads"])

  if(factor <= 0.75 | factor >= 1.5){
    y[, "normalized"] <- round(y[, "reads"]*factor, 3)
    y[y[, "normalized"] < reads.threshold, "normalized"] <- 0
    ratio <- x[, "reads"]/y[, "normalized"]
  } else {
    ratio <- x[, "reads"]/y[, "reads"]
  }
  ratio[is.infinite(ratio) | is.nan(ratio)] <- NA
  ratio.IDs <- paste0(x[, "chr"], ":", x[, "start"])
  ratio.res <- data.frame(x[, c("chr", "start", "end")], ratio)

  if(rm.centromere == TRUE) {
    if(is.null(centromereBins)){
      if(!bin.size %in% c(10000, 25000, 50000, 100000)){
        stop(paste0("SynthEx doesn't have centromere bins pre-calculated for bin size of", bin.size, "; Please use createCentromereBins() to generate the required file or consider to use another bin size."))
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
    ratio.res <- ratio.res[ratio.res[, "chr"]!=23, ]
  }

  res <- list(ratio.res, TRUE)
  names(res) <- c("Ratio", "WGS")
  class(res) <- "WGSRatio"
  return(res)

}
