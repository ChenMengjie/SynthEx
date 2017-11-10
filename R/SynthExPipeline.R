SynthExPipeline <- function(tumor, normal, bin.size, bedTools.dir, genotype.file,
                            result.dir = NULL, working.dir = NULL, verbose = FALSE,
                            saveplot = TRUE, plotNormalized = TRUE,
                            rm.centromere = TRUE, targetAnnotateBins = NULL, centromereBins = NULL, chrX = FALSE,
                            report = TRUE, plot = TRUE, prefix = NULL, reads.threshold = 25, vcf = TRUE,
                            adjust.cutoff = 1.2, seg.count = 200, segmentMethod = "CBS",
                            smoothk = 10, ratio.cutoff = 0.05, K = 1,
                            WGD = 1.35, pos.prop.threhold = 0.6, pos.log2ratio.threhold = 0.75,
                            prop.threshold = 0.0005, delta = 0.1, maf.control = 0.01,
                            tau = 2, sigma = 0.1, len.threshold.K = 10,  group.length.threshold = 2,
                            gain.threshold = log2(1.2), loss.threshold = log2(0.8), lwd = 5){

  ratioCorrectedBias <- SynthExcorrectBias(tumor, normal, bin.size = bin.size, rm.centromere = rm.centromere,
             targetAnnotateBins = targetAnnotateBins, saveplot = saveplot, centromereBins = centromereBins,
             chrX = chrX, plot = plot, result.dir = result.dir, working.dir = working.dir, K =K,
             prefix = prefix, reads.threshold = reads.threshold)

  if(verbose == TRUE) print("Bias correction finished.")

  if(!is.null(genotype.file)){
    ratioNormalized <- normalization(ratioCorrectedBias, bedTools.dir = bedTools.dir, genotype.file = genotype.file, vcf = vcf,
          working.dir = working.dir, result.dir = result.dir, cutoff = reads.threshold, plot = plot, saveplot = saveplot,
          prefix = prefix,  adjust.cutoff = adjust.cutoff, seg.count = seg.count)
    if(verbose == TRUE) print("Normalization finished.")
    ratiotoSeg <- ratioNormalized
  } else {
    ratiotoSeg <- ratioCorrectedBias
  }

  Seg <- createSegments(ratiotoSeg, segmentMethod)

  if(verbose == TRUE) print("Segmentation finished.")

  Segments <- singleCNreport(Seg, report = report, result.dir = result.dir, saveplot = saveplot,
           prefix = prefix, plotNormalized = plotNormalized, WGD = WGD, pos.prop.threhold = pos.prop.threhold,
           pos.log2ratio.threhold = pos.log2ratio.threhold)

  if(!is.null(genotype.file)){
    if(vcf == TRUE){
      PurityCorrected <- purityEstimate(Segments, working.dir = working.dir, result.dir = result.dir, bedTools.dir = bedTools.dir,
      prefix = prefix, report = report, prop.threshold = prop.threshold, delta = delta, maf.control = maf.control,
      tau = tau, sigma = sigma, len.threshold.K = len.threshold.K, group.length.threshold = group.length.threshold,
      gain.threshold = gain.threshold, loss.threshold = loss.threshold, Normalized = plotNormalized)
    } else {
      PurityCorrected <- purityEstimate(Segments, working.dir = working.dir, result.dir = result.dir, bedTools.dir = bedTools.dir,
                                        prefix = prefix, report = report, prop.threshold = prop.threshold, delta = delta, maf.control = maf.control,
                                        tau = tau, sigma = sigma, len.threshold.K = len.threshold.K, group.length.threshold = group.length.threshold,
                                        gain.threshold = gain.threshold, loss.threshold = loss.threshold, Normalized = plotNormalized, vcf = vcf, genotype.file = genotype.file)
    }
    if(verbose == TRUE) print("Purity estimation finished.")
    genomeplot <- chromosomeView(PurityCorrected, prefix = prefix, result.dir = result.dir, saveplot = saveplot, lwd = lwd)
    Segments <- PurityCorrected
  }

  if(verbose == TRUE) print(paste0("Report can be found at: ", result.dir, " ."))

  return(Segments)

}


