.bic.rss <-
  function(RSS, n, Cn, edf, S) 
  {
    ll <- - (log(RSS/n) + log(n)*edf*Cn/n)
    new.r <- ((ll[length(ll)] - ll[-1])/(ll[length(ll)] - 
                                           ll[2])) * (length(ll) - 1) + 1
    diff2 <- diff(new.r, diff = 2) > S
    if (!any(diff2)) 
      return(0)
    maxll = max(which(diff2)) + 1
    return(maxll)
  }

.larsCBS <-
  function(z, selection = .selection.default(), collapse.k = 0, variation.control = TRUE, rss = FALSE, S = 0.1, ...)
  {  
    suppressPackageStartupMessages(require(lars))		
    make.CNA <- function(data){
      chrom <- rep(1, length(data)) 
      maploc <- c(1:length(data)) 
      return(CNA(data, chrom, maploc))
    }
    if(variation.control == TRUE){
      aa <- rnorm(round(length(z)/5), mean=0.5-round(median(z), 1), sd=0.5*sd(z))
      bb <- rnorm(round(length(z)/5), mean=0.5-round(median(z), 1), sd=0.5*sd(z))
      pseudoz <- c(aa, z, bb)		
      y <- pseudoz				
    }
    else{
      y <- z
    }	
    if(collapse.k!=0){
      coly <- collapse(y, collapse.k)	
      oCBS <- segment(make.CNA(coly))
      psi <- c(oCBS$segRows[, 1])[-1]*collapse.k
    }
    else{
      oCBS <- segment(make.CNA(y), verbose = 0, ...)
      psi <- c(oCBS$segRows[, 1])[-1]
    }					
    if (length(psi)==0){ 
      # cat("No estimated breakpoint in the given sequence")	
      result <- list(id.entry = NULL, psiCBS = NULL, psiLARS = NULL, 
                     n.psi = 0, means = round( mean(y), 4), medians = round(median(y), 4))
    }
    else{
      n <- length(y)
      x <- 1:n
      k <- length(psi)
      Z <- matrix(rep(x, k), nrow = n)
      PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
      V <- ifelse((Z > PSI), 1, 0)
      ### Modified from jumpoints() in package cumSeg ### 
      psi0 <- psi
      display1 <- selection$display
      edf.psi <- selection$edf.psi
      type <- selection$type
      Cn <- eval(parse(text = selection$Cn))
      S <- S
      tipoAlg <- selection$alg
      
      olars <- lars(abs(V), y = y, type = tipoAlg, normalize = FALSE, 
                    intercept = TRUE, trace = display1)
      id.var.entry <- (1:ncol(V))[order(olars$entry)]
      if(edf.psi) 
        edf <- (olars$df - 1) * 2 + 1
      else edf <- olars$df  
      
      RSS <- olars$RSS
      if(rss == TRUE)     min.r <- .bic.rss(RSS, n, Cn, edf, S)
      else      min.r <- which.min(log(RSS/n) + log(n) * edf * Cn/n)
      
      id <- sort(c(0, id.var.entry)[1:min.r])  
      id <- id[-1]
      psi1 <- psi[id]
      
      if(variation.control == TRUE){
        psi1 <- psi1[psi1 < (round(length(z)/5)+length(z))-20]-round(length(z)/5)-5
        psi <- psi[psi < (round(length(z)/5)+length(z))-20]-round(length(z)/5)-5	
        psi1 <- psi1[psi1>20]	
        psi <- psi[psi>20]      
        if(collapse.k!=0){
          psi1 <- psi1[psi1 > collapse.k*2 & psi1 < (collapse.k*2+length(z))]
          psi <- psi[psi > collapse.k*2 & psi < (collapse.k*2+length(z))]
        }				
      }	
      start <- c(0, psi1)
      end <- c(psi1, length(z))
      k <- length(start)
      ave <- tapply(1:k, 1:k, function(i){ mean(z[start[i]:end[i]]) })
      med <- tapply(1:k, 1:k, function(i){ median(z[start[i]:end[i]]) })			
      result <- list(id.entry = id.var.entry, psiCBS = psi, psiLARS = psi1, 
                     n.psi = length(psi1), means = round(ave, 4), medians = round(med, 4))
      result$data <- z	          
    }  
    class(result) <- "larsCBSsegmented"
    return(result)
  }

.selection.default <-
  function (display = FALSE, Cn = "log(log(n))", alg = c("stepwise", "lasso"), edf.psi = TRUE) {
    ### Modified from sel.control() in package cumSeg ### 
    alg <- match.arg(alg)
    list(display = display, Cn = Cn, edf.psi = edf.psi, alg = alg)
  }


.tosegment <- function(y){
  diffy <- diff(y)
  changePoints <- c(1:length(y))[diffy != 0]
  changePoints <- changePoints[changePoints!=1 & changePoints != length(y)]
  z <- data.frame(c(1, changePoints), c(changePoints, length(y)), y[c(changePoints, length(y))])
  colnames(z) <- c("start", "end", "value")
  return(z)
}

.chr.tosegment <- function(y, ratioData){
  chrs <- unique(ratioData[, 1])
  allsegments <- NULL
  for(i in 1:length(chrs)){
    suby <- y[ratioData[, 1]%in%chrs[i]]
    subRatio <- ratioData[ratioData[, 1] %in% chrs[i], ]
    segmentSuby <- .tosegment(suby)
    segmentRes <- data.frame(subRatio[segmentSuby[, 1], 1], subRatio[segmentSuby[, 1], 2], subRatio[segmentSuby[, 2], 2], segmentSuby[, 3])
    colnames(segmentRes) <- c("chr", "start", "end", "ratio")
    allsegments <- rbind(allsegments, segmentRes)
  }
  return(allsegments)
}

