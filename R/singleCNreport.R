singleCNreport <- function(Segment, report = FALSE, result.dir = NULL, saveplot = FALSE,
                           prefix = NULL,  plotNormalized = TRUE,
                           WGD = 1.35, pos.prop.threhold = 0.6, pos.log2ratio.threhold = 0.75){

  options(scipen = 50)
  if(is.null(result.dir)) result.dir <- paste0(getwd(), "/result")
  segmentMethod <- Segment$segmentMethod

  if(plotNormalized == TRUE & !is.null(Segment$segmentNormalized)){
    segments <- Segment$segmentNormalized
    data <- Segment$Ratio[, c("chr", "start", "end", "normalized")]
    data$log2ratio <- round(log2(data[, "normalized"] + 0.0001), 4)
  } else {
    segments <- Segment$segmentUnadjusted
    data <- Segment$Ratio[, c("chr", "start", "end", "ratio")]
    data$log2ratio <- round(log2(data[, "ratio"] + 0.0001), 4)
  }

  Segment$WholeGenomeDoubling <- CallWholeGenomeDoubling(segments, WGD = WGD, pos.prop.threhold = pos.prop.threhold, pos.log2ratio.threhold = pos.log2ratio.threhold)$WholeGenomeDoubling

  bin.size <- data[1, "end"] - data[1, "start"] + 1
  maxlength <- tapply(segments$end, segments$chr, max)
  xmax <- sum(maxlength/bin.size) + 400

  if(saveplot == TRUE){
    if(is.null(prefix)){
      jpeg(filename = paste0(result.dir, "/Segmentation_",  segmentMethod, "_", bin.size, ".jpg"), width = 2000, height = 480, quality=100)
    } else {
      jpeg(filename = paste0(result.dir, "/", prefix, "-Segmentation_",  segmentMethod, "_", bin.size, ".jpg"), width = 2000, height = 480, quality=100)
    }
  } else {
    dev.new(width = 20, height = 4.8)
  }
  plot(0, 0, type = "n", ylim = c(-2, 3), xlim = c(0, xmax), axes = F, xlab = "", ylab = "", main = prefix, cex.lab = 2, cex.axis = 2)
  axis(2)
  order <- 0
  for(l in 1:length(unique(segments$chr))){
    sub <- segments[segments$chr == l, ]
    sub.data <- data[data$chr == l, ]
    points(sub.data$start/bin.size+order, sub.data$log2ratio, pch = 20, cex = 0.75)
    for(m in 1:dim(sub)[1]){
      lines(c(sub$start[m]/bin.size+order, sub$end[m]/bin.size+order), c(sub$log2ratio[m], sub$log2ratio[m]), lwd = 5, col = "red")
    }
    order <- order + sub[dim(sub)[1], "end"]/bin.size
    abline(v = order, col = "gray", lty = "dashed", lwd = 2)
    if(l == 23){
      text(order-0.5*sub[dim(sub)[1], "end"]/bin.size, 2.8, "X", font = 2)
    } else {
      text(order-0.5*sub[dim(sub)[1], "end"]/bin.size, 2.8, l, font = 2)
    }
  }
  if(Segment$WholeGenomeDoubling == TRUE){
    abline(h = log2(1/Segment$Adjust), col = "brown", lty = "dashed", lwd = 2)
    text(xmax - 2000, -2, "Whole Genome Doubling", cex = 2, col = "brown")
  }

  if(saveplot == TRUE){
    dev.off()
  }
  if(report == TRUE){
    if(plotNormalized == TRUE & !is.null(Segment$segmentNormalized)){
      segments <- Segment$segmentNormalized
    } else {
      segments <- Segment$segmentUnadjusted
    }

    segments <- segments[, c("chr", "start", "end", "log2ratio")]
    segments <- as.matrix(segments)
    segments[, 1] <- gsub("X", 23, segments[, 1])
    segments <- as.data.frame(segments)
    segments <- segments[segments$chr != "Y", ]

    if(is.null(prefix)){
      dataSample <- segments
      colnames(dataSample) <-c("chr", "start", "end", "seg.mean")
      if(!is.null(prefix)){
        write.table(dataSample, paste0(result.dir, "/", prefix, "-Segment-", segmentMethod, "-", bin.size, ".txt"),
                    col.names = T, row.names = F, sep = "\t", quote = F)
      } else {
        write.table(dataSample, paste0(result.dir, "/Segment-", segmentMethod, "-", bin.size, ".txt"),
                    col.names = T, row.names = F, sep = "\t", quote = F)
      }
    } else {
      dataSample <- cbind(prefix, segments)
      colnames(dataSample) <-c("sample", "chr", "start", "end", "seg.mean")
      if(!is.null(prefix)){
        write.table(dataSample, paste0(result.dir, "/", prefix, "-Segment-", segmentMethod, "-", bin.size, ".txt"),
                  col.names = T, row.names = F, sep = "\t", quote = F)
      } else {
        write.table(dataSample, paste0(result.dir, "/Segment-", segmentMethod, "-", bin.size, ".txt"),
                    col.names = T, row.names = F, sep = "\t", quote = F)
      }
    }
  }
  return(Segment)
}

