create_synthetic_library <- function(counts, foldEnrichment, groupID, interval = 50000000){
  
  Reads <- apply(counts[, -c(1:3)], 2, sum)
  ReadsM <- round(Reads/interval)/(100000000/interval)
  group <- unique(groupID)
  
  syntheticLibrary <- NULL
  syntheticLibrary$Protocol <- list(NULL)
  
  for(i in 1:length(group)){
    ids <- which(groupID == i)
    tab <- table(foldEnrichment[ids], ReadsM[ids])
    seq1 <- as.numeric(colnames(tab))
    seq2 <- as.numeric(rownames(tab))
    syntheticNormal <- create_synthetic_normal(counts[, -c(1:3)][, ids], ReadsM[ids], foldEnrichment[ids], seq1, seq2)
    syntheticLibrary$Protocol[[i]] <- syntheticNormal  
  }
  
  syntheticLibrary$Bins <- counts[, c(1:3)]
  if(substr(syntheticLibrary$Bins[1, 1], 1, 3) == "chr") {
    syntheticLibrary$Bins[, 1] <- gsub("chr", "", syntheticLibrary$Bins[, 1])
  }
  colnames(syntheticLibrary$Bins) <- c("start", "end", "ratio")
  
  syntheticLibrary$NumProtocol <- length(group)
  
  syntheticLibrary$interval <- interval
  
  syntheticLibrary$feature <- c("Reads", "Fold_Enrichment")
  
  class(syntheticLibrary) <- "syntheticLibrary"
  
  return(syntheticLibrary)
}