create_synthetic_normal <- function(counts, normalReads, foldEnrichment, seq1, seq2){
  syntheticNormal <- list(NULL)
  syntheticNormal$readsLevel <- seq1
  syntheticNormal$foldChangeLevel <- seq2
  values <- array(data = NA, dim = c(length(syntheticNormal$readsLevel), length(syntheticNormal$foldChangeLevel), nrow(counts)), dimnames = c("read", "fold", "value"))
  for(i in 1:length(syntheticNormal$readsLevel)){
    for(j in 1:length(syntheticNormal$foldChangeLevel)){
      id <- which(abs(as.numeric(normalReads) - syntheticNormal$readsLevel[i]) < 0.001 & abs(as.numeric(foldEnrichment) - syntheticNormal$foldChangeLevel[j]) < 0.001)
      if(length(id) > 1){
        values[i, j, ] <- apply(counts[, id], 1, median)
      } else if (length(id) == 1){
        values[i, j, ] <- counts[, id]
      }
    }	
  }
  size <- matrix(0, length(syntheticNormal$readsLevel), length(syntheticNormal$foldChangeLevel))
  for(i in 1:length(syntheticNormal$readsLevel)){
    for(j in 1:length(syntheticNormal$foldChangeLevel)){
      id <- which(abs(as.numeric(normalReads) - syntheticNormal$readsLevel[i]) < 0.001 & abs(as.numeric(foldEnrichment) - syntheticNormal$foldChangeLevel[j]) < 0.001)
      size[i, j] <- length(id)
    }
  }
  
  syntheticNormal$values <- values
  syntheticNormal$size <- size
  rownames(syntheticNormal$size) <- paste0("Reads", syntheticNormal$readsLevel)
  colnames(syntheticNormal$size) <- paste0("Fold", syntheticNormal$foldChangeLevel)  
  syntheticNormal$feature <- c("Reads", "Fold_Enrichment")
  class(syntheticNormal) <- "syntheticNormal"
  return(syntheticNormal)
}
