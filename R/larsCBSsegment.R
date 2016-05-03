larsCBSsegment <-
  function(data, selection = .selection.default(), collapse.k = 0, verbose = TRUE, variation.control = FALSE, rss = FALSE, S = 0.1, k = 10,...)
  {
    y <- as.data.frame(data, stringsAsFactors = FALSE)  
    colnames(y) <- c("seqnames", "start", "ratio")
    suppressPackageStartupMessages(require(foreach))
    suppressPackageStartupMessages(require(DNAcopy))
    chrs <- unique(y$seqnames)
    i <- NULL 
    res <- foreach(i=1:length(chrs), .combine=rbind) %do% {		
      sub <- y[y$seqnames==chrs[i], ]
      if(verbose==TRUE){
        message(paste("Processing ", chrs[i], "...", sep=""))
      }
      if(k == 0 | length(sub[, 1]) < 2*k){
        smos <- sub$ratio
      } else {
        smos <- denoise(data = sub$ratio, k = k)
      }
      olarsCBS <- .larsCBS(z = smos, selection = selection, collapse.k = collapse.k, variation.control = variation.control, rss = rss, S = S)# , ...)
      if(olarsCBS$n.psi==0){
        allone <- c(as.character(chrs[i]), sub$start[1], sub$start[length(smos)], olarsCBS$medians)
        names(allone)[1:4] <- c("seqnames", "start", "end", "ratio")
      }
      else{
        one <- cbind(as.character(sub$seqnames[olarsCBS$psiLARS]), sub$start[olarsCBS$psiLARS])
        one <- rbind(c(as.character(sub$seqnames[olarsCBS$psiLARS[1]]), sub$start[1]), one)
        one <- cbind(one, c(one[-1, 2], sub$start[length(smos)]))
        allone <- cbind(one, olarsCBS$medians)
        names(allone)[1:4] <- c("seqnames", "start", "end", "ratio")
      }	
      return(allone)
    }
    rownames(res) <- 1:nrow(res)
    res <- as.data.frame(res)
    return(res)
  }
