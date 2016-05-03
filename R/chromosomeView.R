#' Estimate purity
#'
#' @param Segment
#' @return An object called \code{Segment}

chromosomeView <- function(Segment, prefix = NULL, saveplot = FALSE, result.dir = NULL, lwd = 5){

    if(class(Segment) != "Segment")
    stop("Invalid input class for chromosomeView().")

    segmentMethod <- Segment$sampleinfo$segmentMethod

    if(saveplot == TRUE){
      if(!is.null(prefix)){
        jpeg(paste0(result.dir, "/", prefix, "_purity_corrected-", segmentMethod,  ".jpg"), width = 2000, height = 1600, quality = 100)
      } else {
        jpeg(paste0(result.dir, "/purity_corrected-", segmentMethod,  ".jpg"), width = 2000, height = 1600, quality = 100)
      }
    }

    par(mfrow = c(2, 1), mar = c(2, 2, 2, 2))

    if(Segment$WholeGenomeDoubling == TRUE) {
      ylim <- c(-0.1, 8.5)
    } else {
      ylim <- c(-0.1, 5.5)
    }
    data <- Segment$Ratio[, c("chr", "start", "end")]
    bin.size <- data[1, "end"] -  data[1, "start"] + 1
    segments <- Segment$copynumber[, c("chr", "start", "end", "integerCopy")]
    colnames(segments) <- c("chr", "start", "end", "copy")
    maxlength <- tapply(segments$end, segments$chr, max)
    xmax <- sum(maxlength/bin.size)*1.1

    plot(0, 0, type = "n", ylim = ylim, xlim = c(0, xmax), axes = F, xlab = "", ylab = "", main = paste0("Copy Number \n", prefix), cex.lab = 2, cex.axis = 2)
    axis(2)
    order <- 0
    for(l in c(1:22)){
      sub <- segments[segments$chr == l, ]
      for(m in 1:dim(sub)[1]){
        lines(c(sub$start[m]/bin.size+order, sub$end[m]/bin.size+order), c(sub$copy[m], sub$copy[m]), lwd = lwd)
      }
      order <- order + sub[dim(sub)[1], "end"]/bin.size
      abline(v = order, col = "gray", lty = "dashed", lwd = 2)
      text(order-0.5*sub[dim(sub)[1], "end"]/bin.size, 0.5, l, cex = 0.8)
    }
    if(Segment$WholeGenomeDoubling == TRUE) {
      text(xmax - 2000, ylim[2]-1, "Whole Genome Doubling", cex = 2, col = "brown")
    }

    segments <- Segment$event
    segments$color <- ifelse(segments$event == "Gain", "red",
                      ifelse(segments$event == "Loss", "green", "darkgray"))
    plot(0, 0, type = "n", ylim = c(-2, 3), xlim = c(0, xmax), axes = F, xlab = "", ylab = "", cex.lab = 1.3, cex.axis = 1.3, main = paste0("Event \n", prefix))
    axis(2)
    order <- 0
    for(l in c(1:22)){
      sub <- segments[segments$chr == l, ]
      for(m in 1:dim(sub)[1]){
        lines(c(sub$start[m]/bin.size+order, sub$end[m]/bin.size+order), c(sub$correctedlog2ratiobyPurity[m], sub$correctedlog2ratiobyPurity[m]), lwd = lwd, col = sub$color[m])
      }
      order <- order + sub[dim(sub)[1], "end"]/bin.size
      abline(v = order, col = "gray", lty = "dashed", lwd = 2)
      text(order-0.5*sub[dim(sub)[1], "end"]/bin.size, -1.8, l, cex = 0.8)
    }
    legend("topright", c("Gain", "Loss"), lty = 1, col = c("red", "green"), cex = 2, lwd = lwd)
    if(saveplot == TRUE){
      dev.off()
    }
}



