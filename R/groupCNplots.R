groupCNplots <- function(grps, sample.dirs, gain.thresh = 0.26, loss.thresh = -0.32, species, segmentMethod, returnSummary = F, bin.size = 100000){
  options(scipen = 50)
  totalChrLengths <- read.table(species, sep = "\t")
  threshold <- 0
  segThreshold <- log10(100000)
  cls <- read.table(grps, sep = "\t", na.strings = "", header = T)
  cls.names <- names(table(cls[, 2]))
  samples <- cls[, 1]
  outputFile <- strsplit(grps, ".txt")[1][1]
  data <- NULL
  for(i in 1:dim(samples)[1]){
    sample.name <- samples[i]
    print(paste0("Sample # ", i, ": ", sample.name))
    sampleData <- read.table(paste0(sample.dirs[i], "/SampleSegmentInput-", segmentMethod, ".txt"), as.is = TRUE, sep = "\t", na.strings = "", header = T)
    sampleData <- sampleData[order(sampleData$chr, sampleData$start, sampleData$end), ]
    data <- rbind(data, sampleData)
  }
  
  # Remove segments with logR  >  than gain segment threshold & remove segments with logR < than loss segment threshold
  data <- data[data[, 5]  >  gain.thresh | data[, 5] <  loss.thresh, ]
  alldata <- list()
  
  for(k in 1:length(cls.names)){
    samps <- cls[cls[, 2] == cls.names[k], 1]
    samps <- as.character(samps[!is.na(samps)])
    nsamples <- nlevels(as.factor(samps))
    x <- as.data.frame(data[data[, 1]%in%samps, ])
    if(nrow(x)  >  0){
      vals <- as.numeric(as.vector(t(as.matrix(x)[, 6])))
      size <- as.numeric(as.vector(t(as.matrix(x)[, 4]))) - as.numeric(as.vector(t(as.matrix(x)[, 3])))
      vals[is.na(vals)] <- 0
      size <- log10(size)
      size[is.na(size)] <- 0
      highcut <- threshold
      lowcut <- -1*threshold
      calls <- rep("noChange", nrow(x))
      gains <- rep(-1, nrow(x))
      losses <- rep(-1, nrow(x))
      calls[vals > highcut & size > segThreshold] <- "Gain"
      calls[vals<lowcut & size > segThreshold] <- "Loss"
      gains[vals > highcut & size > segThreshold] <- 1
      losses[vals<lowcut & size > segThreshold] <- 1
      alldata <- cbind(x, calls, gains, losses)
      # format the table and get out some variables
      newTable <- as.matrix(alldata)
      newTable <- as.matrix(apply(newTable, 2, function(x){gsub(" ", "", x)}))
      if(nrow(x) == 1){
        newTable <- as.matrix(t(newTable)) ##if matrix is 1 row, it gets cast to a character vector, which is bad
      }
      samples <- newTable[, 1]
      #nsamples <- nlevels(as.factor(samples))
      #browser()
      if(nrow(x) == 1){
        newTable <- as.numeric(newTable[, -c(1, 7)])
        newTable <- as.matrix(t(newTable)) ##if matrix is 1 row, it gets cast to a character vector, which is bad
      }
      else{
        newTable <- apply(newTable[, c(-1, -7)], 2, as.numeric)
      }
      chr <- newTable[, 1]
      locstop <- newTable[, 3]
      locstart <- newTable[, 2]
      gains <- newTable[, 6]
      losses <- newTable[, 7]
      chr.gains <- vector()
      chr.losses <- vector()
      chr.losses.length <- vector()
      chr.gains.length <- vector()
      chr.losses.start <- vector()
      chr.gains.start <- vector()
      chr.losses.stop <- vector()
      chr.gains.stop <- vector()
      chr.losses.chr <- vector()
      chr.gains.chr <- vector()
      chr.start <- rep(0, 10000)
      chr.end <- rep(0, 10000)
      chr.chr <- rep(0, 10000)
      chr.barsup <- rep(0, 10000)
      chr.barsdown <- rep(0, 10000)
      outable <- matrix(nrow = 10000, ncol = 4)
      gain.first <- 1
      loss.first <- 1
      
      for(tchr in 1:nrow(totalChrLengths)){
        # parse out the data for this chromosome
        this.locstart <- newTable[chr == tchr, 2]
        this.locstop <- newTable[chr == tchr, 3]
        this.gains <- newTable[chr == tchr, 6]
        this.losses <- newTable[chr == tchr, 7]
        this.samples <- samples[chr == tchr]
        
        # sort within sample blocks and merge when contiguous blocks are the same sample and same call
        this.order <- order(this.samples, this.locstart)
        s.gains <- this.gains[this.order]
        s.losses <- this.losses[this.order]
        s.locs.start <- this.locstart[this.order]
        s.locs.stop <- this.locstop[this.order]
        s.samples <- this.samples[this.order]
        indices.start <- rep(0, 10000)
        indices.stop <- rep(0, 10000)
        count <- 1
        index <- 1
        if(length(s.locs.start) > 1){
          while((length(s.locs.start)-1)  >  index){
            tstart <- index
            looped <- F
            while(s.samples[index] == s.samples[index+1] && s.gains[index] == s.gains[index+1] && s.losses[index] == s.losses[index+1] && (length(s.locs.start)-1)  >  index){
              index <- index+1
              looped <- T
            }
            if(looped){index <- index-1}
            indices.start[count] <- tstart
            indices.stop[count] <- index
            count <- count+1
            index <- index+1
          }
          if(sum(indices.start > 0) > 0){
            s.gains <- s.gains[indices.start]
            s.losses <- s.losses[indices.start]
            s.locs.start <- s.locs.start[indices.start]
            s.locs.stop <- s.locs.stop[indices.stop]
          }
        }
        
        this.order <- order(s.locs.start)
        
        # sort across locations
        s.gains <- s.gains[this.order]
        s.losses <- s.losses[this.order]
        s.locs.start <- s.locs.start[this.order]
        s.locs.stop <- s.locs.stop[this.order]
        s.samples <- s.samples[this.order]
        
        # find the overlaping gain segments
        gain.filter <- s.gains == 1
        if(sum(gain.filter) > 0){
          sg.locs.start <- s.locs.start[gain.filter]
          sg.locs.stop <- s.locs.stop[gain.filter]
          
          breakType <- c(rep("start", length(sg.locs.start)), rep("stop", length(sg.locs.stop)))
          breakLocs <- c(sg.locs.start, sg.locs.stop)
          breakType <- breakType[sort.list(breakLocs)]
          breakLocs <- breakLocs[sort.list(breakLocs)]		
          
          segCount <- 1
          segStop <- vector()
          segStart <- vector()
          segScores <- vector()
          segStart[1] <- breakLocs[1]
          segScore <- 1
          lastStart <- segStart[1]
          lastStop <- 0
          for(j in 2:length(breakLocs)){
            if(breakType[j] == "start"){
              if(breakLocs[j]  ==  lastStart){
                segScore <- segScore + 1
              }else{
                if(segScore > 0){
                  segStop[segCount] <- breakLocs[j]
                  lastStop <- segStop[segCount]
                  segScores[segCount] <- segScore
                  segCount <- segCount + 1
                }
                segStart[segCount] <- breakLocs[j]
                segScore <- segScore + 1
                lastStart <- segStart[segCount]
              }
              
            }else{
              if(breakLocs[j]  ==  lastStop){
                segScore <- segScore - 1
              }else{
                segStop[segCount] <- breakLocs[j]
                lastStop <- segStop[segCount]
                segScores[segCount] <- segScore
                segCount <- segCount + 1 	
                segScore <- segScore - 1
                segStart[segCount] <- breakLocs[j]
              }
            }
          }
          gain.start <- segStart[-1*length(segStart)]
          gain.stop <- segStop
          gain <- segScores
        }else{
          gain.start <- 0
          gain.stop <- totalChrLengths[totalChrLengths[, 1] == tchr, 2]
          gain <- 0
        }		
        
        # find the overlaping loss segments
        loss.filter <- s.losses == 1
        if(sum(loss.filter) > 0){
          sl.locs.start <- s.locs.start[loss.filter]
          sl.locs.stop <- s.locs.stop[loss.filter]
          
          breakType <- c(rep("start", length(sl.locs.start)), rep("stop", length(sl.locs.stop)))
          breakLocs <- c(sl.locs.start, sl.locs.stop)
          breakType <- breakType[sort.list(breakLocs)]
          breakLocs <- breakLocs[sort.list(breakLocs)]		
          
          segCount <- 1
          segStop <- vector()
          segStart <- vector()
          segScores <- vector()
          segStart[1] <- breakLocs[1]
          segScore <- 1
          lastStart <- segStart[1]
          lastStop <- 0
          for(j in 2:length(breakLocs)){
            if(breakType[j] == "start"){
              if(breakLocs[j]  ==  lastStart){
                segScore <- segScore + 1
              }else{
                if(segScore > 0){
                  segStop[segCount] <- breakLocs[j]
                  lastStop <- segStop[segCount]
                  segScores[segCount] <- segScore
                  segCount <- segCount + 1
                }
                segStart[segCount] <- breakLocs[j]
                segScore <- segScore + 1
                lastStart <- segStart[segCount]
              }
              
            }else{
              if(breakLocs[j]  ==  lastStop){
                segScore <- segScore - 1
              }else{
                segStop[segCount] <- breakLocs[j]
                lastStop <- segStop[segCount]
                segScores[segCount] <- segScore
                segCount <- segCount + 1 	
                segScore <- segScore - 1
                segStart[segCount] <- breakLocs[j]
              }
            }
          }
          loss.start <- segStart[-1*length(segStart)]
          loss.stop <- segStop
          loss <- segScores
        }else{
          loss.start <- 0
          loss.stop <- totalChrLengths[totalChrLengths[, 1] == tchr, 2]
          loss <- 0
        }		
        
        # finalize the segments into the entire chrom
        loss.index <- loss.first
        if(loss.start[1] !=  0){
          chr.losses[loss.index] <- 0
          chr.losses.length[loss.index] <- loss.start[1]
          chr.losses.start[loss.index] <- 0
          chr.losses.stop[loss.index] <- loss.start[1]
          loss.index <- loss.index+1
        }
        chr.losses[loss.index] <- loss[1]
        chr.losses.length[loss.index] <- loss.stop[1]-loss.start[1]
        chr.losses.start[loss.index] <- loss.start[1]
        chr.losses.stop[loss.index] <- loss.stop[1]
        loss.index <- loss.index+1
        if(length(loss.start) > 1){
          for(j in 2:length(loss.start)){
            if(loss.start[j]!= loss.stop[j-1]){
              chr.losses[loss.index] <- 0
              chr.losses.length[loss.index] <- loss.start[j] - loss.stop[j-1]
              chr.losses.start[loss.index] <- loss.stop[j-1]
              chr.losses.stop[loss.index] <- loss.start[j]
              loss.index <- loss.index + 1
            }
            chr.losses[loss.index] <- loss[j]
            chr.losses.length[loss.index] <- loss.stop[j]-loss.start[j]
            chr.losses.start[loss.index] <- loss.start[j]
            chr.losses.stop[loss.index] <- loss.stop[j]
            loss.index <- loss.index + 1
          }
        }
        if(totalChrLengths[tchr, 2] > sum(chr.losses.length[loss.first:(loss.index-1)])){
          chr.losses[loss.index] <- 0
          chr.losses.start[loss.index] <- sum(chr.losses.length[loss.first:(loss.index-1)])
          chr.losses.stop[loss.index] <- totalChrLengths[totalChrLengths[, 1] == tchr, 2]
          chr.losses.length[loss.index] <- totalChrLengths[totalChrLengths[, 1] == tchr, 2]-sum(chr.losses.length[loss.first:(loss.index-1)])
          loss.index <- loss.index + 1
        }
        if(totalChrLengths[tchr, 2]<sum(chr.losses.length[loss.first:(loss.index-1)])){
          print(paste("Warning: Lenghts are longer than total chrom length for chrom", tchr))
          print(paste("					Chrom is", totalChrLengths[tchr, 2], "bases"))
          print(paste("					Calculated length is", sum(chr.losses.length[loss.first:(loss.index-1)])))
          
        }
        chr.losses.chr[loss.first:(loss.index-1)] <- rep(tchr, loss.index-loss.first)
        loss.first <- loss.index		
        
        gain.index <- gain.first
        if(gain.start[1] !=  0){
          chr.gains[gain.index] <- 0
          chr.gains.length[gain.index] <- gain.start[1]
          chr.gains.start[gain.index] <- 0
          chr.gains.stop[gain.index] <- gain.start[1]
          gain.index <- gain.index+1
        }
        chr.gains[gain.index] <- gain[1]
        chr.gains.length[gain.index] <- gain.stop[1]-gain.start[1]
        chr.gains.start[gain.index] <- gain.start[1]
        chr.gains.stop[gain.index] <- gain.stop[1]
        gain.index <- gain.index+1
        if(length(gain.start) > 1){
          for(j in 2:length(gain.start)){
            if(gain.start[j]!= gain.stop[j-1]){
              chr.gains[gain.index] <- 0
              chr.gains.length[gain.index] <- gain.start[j] - gain.stop[j-1]
              chr.gains.start[gain.index] <- gain.stop[j-1]
              chr.gains.stop[gain.index] <- gain.start[j]
              gain.index <- gain.index + 1
            }
            chr.gains[gain.index] <- gain[j]
            chr.gains.length[gain.index] <- gain.stop[j]-gain.start[j]
            chr.gains.start[gain.index] <- gain.start[j]
            chr.gains.stop[gain.index] <- gain.stop[j]
            gain.index <- gain.index + 1
          }
        }
        if(totalChrLengths[totalChrLengths[, 1] == tchr, 2] > sum(chr.gains.length[gain.first:(gain.index-1)])){
          chr.gains[gain.index] <- 0
          chr.gains.start[gain.index] <- sum(chr.gains.length[gain.first:(gain.index-1)])
          chr.gains.stop[gain.index] <- totalChrLengths[totalChrLengths[, 1] == tchr, 2]
          chr.gains.length[gain.index] <- totalChrLengths[totalChrLengths[, 1] == tchr, 2]-sum(chr.gains.length[gain.first:(gain.index-1)])
          gain.index <- gain.index + 1
        }
        if(totalChrLengths[tchr, 2] < sum(chr.gains.length[gain.first:(gain.index-1)])){
          print(paste("Warning: Lenghts are longer than total chrom length for chrom", tchr))
          print(paste("					Chrom is", totalChrLengths[tchr, 2], "bases"))
          print(paste("					Calculated length is", sum(chr.gains.length[gain.first:(gain.index-1)])))
          
        }
        chr.gains.chr[gain.first:(gain.index-1)] <- rep(tchr, gain.index-gain.first)
        gain.first <- gain.index
      }
      
      jpeg(filename = paste(cls.names[k], "_SWITCHdna_landscape.jpg", sep = ""), width = 2000, height = 480, quality = 100)
      barplot((chr.gains/nsamples), chr.gains.length, space = 0, col = "darkred", border = NA, ylim = c(-1, 1), main = cls.names[k])
      barplot((-1*chr.losses/nsamples), chr.losses.length, space = 0, col = "darkgreen", border = NA, add = T)
      
      count <- 1
      chr.lengths <- vector()
      chr.pos = 0
      for(j in 1:nrow(totalChrLengths)){
        if(j > 1){
          chr.pos = sum(as.numeric(totalChrLengths[1:(j-1), 2]))
        }
        else{
          chr.pos = 0
        }
        chr.lengths[count] <- totalChrLengths[j, 2]-50
        
        chr.barsup[count] <- 0
        chr.barsdown[count] <- 0
        count <- count+1
        chr.lengths[count] <- 50
        chr.barsup[count] <- 1
        chr.barsdown[count] <- -1			
        count <- count+1
        abline(v = chr.pos, col = "gray", lwd = 3)
        text(chr.pos, -0.9, labels = j, col = "black", cex = 1.5)
      }
      
      dev.off()
      
      outable <- rbind(cbind(chr.gains.chr, chr.gains.start, chr.gains.stop, (chr.gains/nsamples)), cbind(chr.losses.chr, chr.losses.start, chr.losses.stop, (-1*chr.losses/nsamples)))
      outable <- as.matrix(outable[!(is.na(outable[, 4]) | outable[, 4] == 0), ])
      if(ncol(outable) == 1){
        outable <- t(outable) ##this exists because if there is only one region for a summary, then it gets cast to a vector and not the 4-column matrix it is supposed to
      }
      if(length(outable) > 0){
        outable <- rbind(c("Chr", "Start", "End", "Percent Change"), outable)
        write.table(outable, paste(cls.names[k], "_SWITCHdna_summary.txt", sep = ""), sep = "\t", col.names = F, row.names = F)
      }
      else{
        print(paste(cls.names[k], "_SWITCHdna_summary.txt contains no elements!", sep = ""))
      }
    }
    else{
      print(paste(cls.names[k], "_SWITCHdna_summary.txt contains no elements after filter!", sep = ""))
    }
  }
  if(returnSummary == T){
    return(outable)
  }
  
}
