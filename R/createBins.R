createBins <- function(binsize = 100000, chromLength = NULL, write = TRUE, result.dir = NULL, prefix = "hg19"){
  options(scipen = 50)
  if(!is.null(chromLength)){
    x <- read.table(chromLength, as.is = T)
  } else {
    load(hg19.chromosome.size)
    x <- hg19.chromosome.size
  }
  chromosomes <- x[, 2]
  allintervals <- NULL
  for(i in c(1:24)){
    max <- ceiling(chromosomes[i]/binsize)
    start <- seq(1, by = binsize, length.out = max)
    end <- start + binsize - 1
    res <- cbind(x[i, 1], start, end)
    allintervals <- rbind(allintervals, res)
  }
  if(write == TRUE){
    write.table(allintervals, paste0(result.dir, "/", prefix, "_", binsize, ".bed"), col.names = F, row.names = F, quote = F, sep = "\t")
  }
  return(allintervals)
}


