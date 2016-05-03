createCentromereBins <- function(cytoband = NULL, bin.size = 100000, write = FALSE, result.dir = NULL){
  options(scipen = 50)
  if(is.null(cytoband)){
    data(CentromereAnnotations)
    cytoband <- CentromereAnnotations$cytoband
  }

  cytoband[, 1] <- gsub("chr", "", cytoband[, 1])
  centromere <- NULL
  for(i in c(1:22, "X")){
    sub.cytoband <- cytoband[cytoband[, 1] == i, ]
    pq <- substr(sub.cytoband[, 4], 1, 1)
    pq.neighbor <- paste0(pq[-length(pq)], pq[-1])
    window <- which(pq.neighbor == "pq")
    centromere <- rbind(centromere, sub.cytoband[window:(window + 1), ])
  }
  allbins <- NULL
  for(j in 1:nrow(centromere)){
    start <- seq(centromere[j, 2], centromere[j, 3] - bin.size, bin.size)
    allbins <- rbind(allbins, cbind(centromere[j, 1], start + 1, start + bin.size))
  }
  if(write == TRUE){
    write.table(allbins, paste0(result.dir, "/centromere.bin", binsize, ".bed"), col.names = F, row.names = F, quote = F, sep ="\t")
  }
  allbins <- as.data.frame(allbins)
  return(allbins)
}






