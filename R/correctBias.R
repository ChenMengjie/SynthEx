correctBias <- function(tumor, normal, bin.size = 100000, rm.centromere = TRUE,
                        targetAnnotateBins = NULL, centromereBins = NULL, chrX = FALSE, plot = TRUE,
                        saveplot = TRUE, result.dir = NULL, prefix = NULL, reads.threshold = 50){

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
    centromere[, 2] <- centromere[, 2] + 1
    centromere.IDs <- paste0(centromere[, 1], ":", centromere[, 2])
    ratio.IDs <- paste0(ratio.res[, "chr"], ":", ratio.res[, "start"])
    ratio.res <- ratio.res[! ratio.IDs%in%centromere.IDs, ]
  }

  ratio <- ratio.res[, "ratio"]
  ratio.IDs <- paste0(ratio.res[, "chr"], ":", ratio.res[, "start"])
  if(is.null(targetAnnotateBins)){
    if(!bin.size %in% c(10000, 25000, 50000, 100000)){
      stop(paste0("SynthEx doesn't have centromere bins pre-calculated for bin size of", bin.size, "; Please use createTargetBins()
                  to generate the required file or consider to use another bin size."))
    } else {
      data(TargetAnnotations)
      ss <- paste0("target <- TargetAnnotations$bin", bin.size)
      eval(parse(text=ss))
    }
  } else {
    target <- read.delim(targetAnnotateBins, header = F, as.is = T)
  }

  target[, 5] <- target[, 5]+1
  target.IDs <- paste0(target[, 4], ":", target[, 5])
  target.IDs <- target.IDs[!duplicated(target.IDs)]

  #identify target and off target regions
  target.status <- ifelse(ratio.IDs %in% target.IDs, "target", "non-target")

  all.dis.target.left <- NULL
  all.dis.target.right <- NULL
  all.dis.non.target.left <- NULL
  all.dis.non.target.right <- NULL
  all.dis.within.target.left <- NULL
  all.dis.within.target.right <- NULL

  allchrs <- ifelse(chrX == TRUE, 1:23, 1:22)
  for(i in allchrs){
    chr.ratio <- ratio[grep(paste0(i, ":"), ratio.IDs)]
    chr.names <- ratio.IDs[grep(paste0(i, ":"), ratio.IDs)]
    chr.pos <- unlist(strsplit(chr.names, ":"))[seq(2, 2*length(chr.names), by=2)]
    chr.pos <- (as.numeric(chr.pos)-1)/bin.size
    chr.target.status <- target.status[grep(paste0(i, ":"), ratio.IDs)]
    chr.ratio <- chr.ratio[order(chr.pos)]
    target.bins <- c(1:length(chr.ratio))[chr.target.status == "target"]
    all.neighbor.bins <- sort(unique(c(target.bins - 1, target.bins, target.bins + 1)))
    sub.chr.ratio <- chr.ratio[all.neighbor.bins]

    #Calculate Target-NonTarget Ratios
    select.target.bins <- target.bins[!target.bins %in% c(target.bins - 1) & !target.bins %in% c(target.bins + 1, 1, nrow(chr.ratio))]
    select.target.bins.ratio <- chr.ratio[select.target.bins]
    select.left.target.bins.ratio <- chr.ratio[select.target.bins - 1]
    select.right.target.bins.ratio <- chr.ratio[select.target.bins + 1]
    dis.target.left <- select.target.bins.ratio - select.left.target.bins.ratio
    dis.target.right <- select.target.bins.ratio - select.right.target.bins.ratio

    #Calculate NonTarget-NonTarget ratios
    non.target.bins <- c(1:length(chr.ratio))[-c(target.bins, 1, nrow(chr.ratio))]
    non.target.bins.ratio <- chr.ratio[non.target.bins]
    non.left.target.bins.ratio <- chr.ratio[non.target.bins - 1]
    non.right.target.bins.ratio <- chr.ratio[non.target.bins + 1]
    dis.non.target.left <- non.target.bins.ratio - non.left.target.bins.ratio
    dis.non.target.right <- non.target.bins.ratio- non.right.target.bins.ratio

    #Calculate Target-Target Ratios
    select.left.target.target.bins <- target.bins[target.bins %in% c(target.bins - 1)]
    select.right.target.target.bins <- target.bins[target.bins %in% c(target.bins + 1)]
    dis.target.target.left <- chr.ratio[select.left.target.target.bins] - chr.ratio[select.left.target.target.bins + 1]
    dis.target.target.right <- chr.ratio[select.right.target.target.bins] - chr.ratio[select.right.target.target.bins - 1]

    all.dis.target.left <- c(all.dis.target.left, dis.target.left)
    all.dis.target.right <- c(all.dis.target.right, dis.target.right)
    all.dis.non.target.left <- c(all.dis.non.target.left, dis.non.target.left)
    all.dis.non.target.right <- c(all.dis.non.target.right, dis.non.target.right)
    all.dis.within.target.left <- c(all.dis.within.target.left, dis.target.target.left)
    all.dis.within.target.right <- c(all.dis.within.target.right, dis.target.target.right)

  }

  all.dis.target <- c(all.dis.target.left, all.dis.target.right)
  all.dis.non.target <- c(all.dis.non.target.left, all.dis.non.target.right)
  all.dis.within.target <- c(all.dis.within.target.left, all.dis.within.target.right)

  d1 <- density(all.dis.target, na.rm = T)
  d2 <- density(all.dis.non.target, na.rm = T)
  d3 <- density(all.dis.within.target.left, na.rm = T)

  all.dis.others <- c(all.dis.non.target, all.dis.within.target)

  target.data <- all.dis.target[!is.nan(all.dis.target) & !is.infinite(all.dis.target)]
  nontarget.data <- all.dis.others[!is.nan(all.dis.others) & !is.infinite(all.dis.others)]

  target.sd <- sd(target.data, na.rm = T)
  non.target.sd <- sd(nontarget.data, na.rm = T)
  target.max <- mean(target.data, na.rm = T)
  target.bias.statistics <- round(c(target.max, mean(nontarget.data, na.rm = T), target.sd, non.target.sd), 3)
  names(target.bias.statistics) <- c("TargetMean", "NonTargetMean", "TargetSd", "NonTargetSd")

  if(plot == TRUE) {
    if(abs(target.max) > 1){
      print("Caution: The ratios between bins at the boundaries of target and non-target regions deviate a lot from zero, indicating the ratios from the target and non-target regions are very different")
      if(saveplot == TRUE){
        if(!is.null(prefix)){
          pdf(paste0(result.dir, "/",  prefix, "-DensityPlot--CautionSample.pdf"))
        } else {
          pdf(paste0(result.dir, "/DensityPlot--CautionSample.pdf"))
        }
      }
      if(is.null(prefix)){
        sample.name <- ""
      } else {
        sample.name <- prefix
      }
      plot(d1, col = "red", xlim = c(-3, 3), ylim = c(0, max(d1$y, d2$y, d3$y)), main = prefix, xlab = "Successive Differenc", lwd = 2)
      lines(d2, lwd = 2)
      lines(d3, col = "blue", lwd = 2)
      legend("topright", c("Target-Non", "Non-Non", "Target-Target"), lty = 1, lwd = 2, col = c("red", "black", "blue"))
      if(saveplot == TRUE){
        dev.off()
      }
    } else{
      if(saveplot == TRUE){
        if(!is.null(prefix)){
          pdf(paste0(result.dir, "/",  prefix, "-DensityPlot.pdf"))
        } else {
          pdf(paste0(result.dir, "/DensityPlot.pdf"))
        }
      }
      if(is.null(prefix)){
        sample.name <- ""
      } else {
        sample.name <- prefix
      }
      plot(d1, col = "red", xlim = c(-3, 3), ylim = c(0, max(d1$y, d2$y, d3$y)), main = sample.name, xlab = "Successive Difference", lwd = 2)
      lines(d2, lwd = 2)
      lines(d3, col = "blue", lwd = 2)
      legend("topright", c("Target-Non", "Non-Non", "Target-Target"), lty = 1, lwd = 2, col = c("red", "black", "blue"))
      if(saveplot == TRUE){
        dev.off()
      }
    }
  }

  #Adjustment factor
  bias <- target.max
  ratio <- ifelse(!is.na(ratio) & !is.nan(ratio) & ratio.IDs %in% target.IDs, ratio - bias, ratio)
  ratio <- ifelse(!is.na(ratio) & !is.nan(ratio) & ratio < 0, 0, ratio)
  ratio[is.infinite(ratio) | is.nan(ratio)] <- NA
  ratio <- round(ratio/median(ratio, na.rm = T), 3)
  ratio.res[, "ratio"] <- ratio

  if(chrX == FALSE){
    ratio.res <- ratio.res[ratio.res[, "chr"]!=23, ]
  }

  res <- list(target.bias.statistics, ratio.res, FALSE)
  names(res) <- c("TargetBiasStatistics", "Ratio", "Synthetic")
  class(res) <- "RatioCorrectBiasInTargets"
  return(res)
}

