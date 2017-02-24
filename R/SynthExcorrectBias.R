SynthExcorrectBias <- function(tumor, normal, bin.size = 100000, rm.centromere = TRUE, K = 1,
                             targetAnnotateBins = NULL, saveplot = FALSE, centromereBins = NULL, chrX = FALSE,
                             plot = TRUE, result.dir = NULL, working.dir = NULL, prefix = NULL, reads.threshold = 25){
  options(scipen = 50)
  if(is.null(result.dir)) result.dir <- "result"
  if(is.null(working.dir)) working.dir <- "working"
  system(paste0("mkdir ", result.dir))
  system(paste0("mkdir ", working.dir))

  if(class(tumor) == "character"){
    tumor <- read.delim(tumor, header = F, as.is = T)
    if(ncol(tumor) != 4)
    stop("Wrong input for tumor sample. A bed file is expected")
  }
  colnames(tumor) <- c("chr", "start", "end", "reads")
  if(substr(tumor[1, 1], 1, 3) == "chr") {
    tumor[, 1] <- gsub("chr", "", tumor[, 1])
  }
  tumor[, 1] <- gsub("X", "23", tumor[, 1])
  tumor[, 1] <- gsub("Y", "24", tumor[, 1])
  len <- tumor[1, "end"] - tumor[1, "start"]
  if(len != bin.size ){
    stop("\"bin.size\" should match with the input file!")
  }
  tumor[, 2] <- tumor[, 2] + 1


  if(class(normal) == "character"){
    normal <- read.delim(normal, header = F, as.is = T)
  }
  if(class(normal) == "data.frame"){
    if(ncol(normal) == 4){
      if(substr(normal[1, 1], 1, 3) == "chr") {
        normal[, 1] <- gsub("chr", "", normal[, 1])
      }
      normal[, 1] <- gsub("X", "23", normal[, 1])
      normal[, 1] <- gsub("Y", "24", normal[, 1])
      colnames(normal) <- c("chr", "start", "end", "reads")
      len2 <- normal[1, "end"] - normal[1, "start"]
      if( len2 != bin.size | len != len2){
        stop("Wrong input for normal sample. Normal can be either one matched sample or any number of multiple samples.")
      }
    }

    normal[, 2] <- normal[, 2] + 1

    if(ncol(normal) > 4){
      if(K > 1){
        Corrected <- synthetic_correctBias_nearsamples(tumor, normal[, -c(1:3)], bin.size = bin.size, rm.centromere = rm.centromere, K = K,
                                                   targetAnnotateBins = targetAnnotateBins, saveplot = saveplot, centromereBins = centromereBins, chrX = chrX,
                                                    plot = plot, result.dir = result.dir, prefix = prefix, reads.threshold = reads.threshold)
      } else {
        Corrected <- synthetic_correctBias_allsamples(tumor, normal[, -c(1:3)], bin.size = bin.size, rm.centromere = rm.centromere,
                                                  targetAnnotateBins = targetAnnotateBins, saveplot = saveplot, centromereBins = centromereBins, chrX = chrX,
                                                  plot = plot, result.dir = result.dir, prefix = prefix, reads.threshold = reads.threshold)
      }
    } else {
      Corrected <- correctBias(tumor, normal, bin.size = bin.size, rm.centromere = rm.centromere,
                             targetAnnotateBins = targetAnnotateBins, saveplot = saveplot, centromereBins = centromereBins, chrX = chrX,
                             plot = plot, result.dir = result.dir, prefix = prefix, reads.threshold = reads.threshold)

    }
  } else if(class(normal) == "syntheticLibrary"){
    Corrected <- synthetic_correctBias_allpossibility(tumor, normal, bin.size = bin.size, rm.centromere = rm.centromere,
                             targetAnnotateBins = targetAnnotateBins, saveplot = saveplot, centromereBins = centromereBins, chrX = chrX,
                             plot = plot, result.dir = result.dir, prefix = prefix, reads.threshold = reads.threshold)
  }

  return(Corrected)
}
