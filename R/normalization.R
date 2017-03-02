normalization <- function(ratioCorrectedBias, bedTools.dir, genotype.file, vcf = TRUE, working.dir = NULL,
                          result.dir = NULL, cutoff = 20, plot = TRUE, saveplot = FALSE,
                          prefix = NULL, adjust.cutoff = 1.2, seg.count = 200){

  freebayes_to_bed <- system.file("doc", "freebayes_to_bed.py", package="SynthEx")

  options(scipen = 50)
  if(is.null(result.dir)) result.dir <- "result"
  if(is.null(working.dir)) working.dir <- "working"

  if(class(ratioCorrectedBias) != "RatioCorrectBiasInTargets" & class(ratioCorrectedBias) != "WGSRatio" ){
    stop("The input of normalization() must be the output from correctBias() or calratioWGS().")
  }

  write.table(ratioCorrectedBias$Ratio, paste0(working.dir, "/Ratio.bed"),
              col.names = F, row.names = F, sep = "\t", quote = F)

  if(vcf==TRUE) {
    ######### calculate MAF from VCF file #########
    ff <- paste0("cat ", genotype.file, " | python ", freebayes_to_bed, " ", cutoff, " 0.05  > ", working.dir, "/tumor.MAF.highcut.bed")
    system(ff)
    ######## Creates the intersect(using intersectBed) of MAF + adjusted bin ratios
    ff <- paste0(bedTools.dir, " -a ", working.dir, "/tumor.MAF.highcut.bed -b ", paste0(working.dir, "/Ratio.bed"),
               " -wa -wb > ",  working.dir, "/tumor.MAF.ratio.bed")
    system(ff)
  } else {
    ff <- paste0(bedTools.dir, " -a ", genotype.file, " -b ", paste0(working.dir, "/Ratio.bed"),
                 " -wa -wb > ",  working.dir, "/tumor.MAF.ratio.bed")
    system(ff)
  }
  ##########
  data <- read.delim(paste0(working.dir, "/tumor.MAF.ratio.bed"), header = F)
  if(substr(data[1, 1], 1, 3) == "chr") {
    data[, 1] <- gsub("chr", "", data[, 1])
  }
  data <- data[data[, 1] %in% c(1:22), ]
  segments <- paste0(data[, 5], ":", data[, 6])
  median.MAF <- tapply(data[, 4], segments, median)
  median.MAF.segments <- data.frame(names(median.MAF), median.MAF)
  colnames(median.MAF.segments) <- c("segments", "MAF")
  ratio.segments <- data.frame(segments, data[, 8])
  ratio.segments <- ratio.segments[!duplicated(ratio.segments), ]
  colnames(ratio.segments) <- c("segments", "ratio")
  merged.segments <- merge(median.MAF.segments, ratio.segments, by.x = "segments", by.y = "segments")
  merged.segments <- merged.segments[!is.na(merged.segments[, 2]) & !is.na(merged.segments[, 3]) & is.finite(merged.segments[, 3]), ]

  all.median.ratio <- merged.segments[merged.segments$MAF > 0.45, "ratio"]
  all.median.AF <- merged.segments[merged.segments$MAF > 0.45, "MAF"]

  adjust <- median(all.median.ratio, na.rm = T)

  median.ratio <- all.median.ratio[all.median.ratio < 2.5*adjust]
  median.AF <- all.median.AF[all.median.ratio < 2.5*adjust]

  rr <- density(median.ratio, na.rm = T)
  adjust2 <- rr$x[which.max(rr$y)]
  sorted.median.ratio <- sort(median.ratio)

  if( ! all(sorted.median.ratio == 1) ){

    suppressPackageStartupMessages(require(mclust))
    #Identify 2n cluster - Diploid reference#
    mccluster <- Mclust(sorted.median.ratio, G = 1:5, prior = priorControl(functionName = "defaultPrior", shrinkage = 0.1))
    means <- mccluster$parameters$mean
    prop <- mccluster$parameters$pro

    prop <- prop[order(means)]
    means <- sort(means)

    #Adjustment

    adjust3 <- NA
    if (length(means) == 2){
      adjust3 <- ifelse(prop[1] > 0.1 & abs(means[1] - adjust2) >= 0.2, means[1], NA)
    } else if (length(means) > 2){
     adjust3 <- ifelse((prop[1] > 0.1 & abs(means[1] - adjust2) >= 0.2) | (prop[1] > 0.5), means[1], NA)
    }
    if(length(means) > 2 & prop[1] < 0.1 & means[2] < 1){
      adjust3 <- ifelse(abs(means[2] - means[3]) >= 0.2, means[2], NA)
    }
    if(length(median.ratio[median.ratio == 0]) > 50){
      adjust2 <- 1
    } else if(!is.na(adjust3)) {
      adjust2 <- adjust3
    }

    if(adjust2 > adjust.cutoff | adjust2 > adjust) adjust2 <- adjust

    if(length(median.ratio) < seg.count) adjust2 <- 1

    adjust2 <- as.numeric(round(adjust2, 3))
  } else {
    adjust2 <- 1
  }

  AF.ratio <- data.frame(median.AF, median.ratio, row.names = NULL)

  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(gridExtra))
  theme_set(theme_bw())

  if(plot == TRUE) {
    #scatterplot of x and y variables
    if(saveplot == TRUE){
      if(!is.null(prefix)){
        jpeg(filename = paste0(result.dir, "/", prefix, "-BaselinePlot.jpg"), width = 1200, height = 600, quality=100)
      } else {
        jpeg(filename = paste0(result.dir, "/BaselinePlot.jpg"), width = 1200, height = 600, quality=100)
      }
    }
    scatter <- ggplot(AF.ratio, aes(median.AF, median.ratio)) +  geom_point() +  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) + xlab("median MAF") + ylab("Tumor/Normal Ratio")
    ypos <- max(density(median.ratio)$y)*0.6
    adjust.data <- data.frame(adjust2, ypos)
    #marginal density of y - plot on the right
    plot_right <- ggplot(AF.ratio, aes(median.ratio, fill = "red")) + geom_density(alpha = .5) + coord_flip() + theme(legend.position = "none")  + geom_vline(xintercept = adjust2, color = "blue", size = 2) + ylab("Density") + xlab("Tumor/Normal Ratio") + geom_text(aes(adjust2 + 0.3, ypos, label = "Baseline: \n Copy Number = 2"), data = adjust.data)
    grid.arrange(scatter, plot_right, ncol = 2, nrow = 1, widths = c(2, 2))
    if(saveplot == TRUE){
      dev.off()
    }
  }

  ratioCorrectedBias$Ratio[, "normalized"] <- round(ratioCorrectedBias$Ratio[, "ratio"]/adjust2, 3)
  if(class(ratioCorrectedBias) == "RatioCorrectBiasInTargets"){
    ratioNormalized <- list(ratioCorrectedBias$TargetBiasStatistics, ratioCorrectedBias$Ratio, adjust2, ratioCorrectedBias$Synthetic)
    names(ratioNormalized) <- c("TargetBiasStatistics", "Ratio", "Adjust", "Synthetic")
  } else {
    ratioNormalized <- list(ratioCorrectedBias$Ratio, adjust2, ratioCorrectedBias$Synthetic, TRUE)
    names(ratioNormalized) <- c("Ratio", "Adjust", "Synthetic", "WGS")
  }

  class(ratioNormalized) <- "RatioNormalized"
  return(ratioNormalized)

}
