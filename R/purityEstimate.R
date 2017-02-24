#' Estimate purity
#'
#' @param Segment
#' @param result.dir The directory to save the results.
#' @param bedTools.dir
#' @param prefix
#' @param prop.threshold
#' @param delta
#' @param maf.control
#' @param tau
#' @param sigma
#' @return An object called \code{Segment}

purityEstimate <- function(Segment, working.dir = NULL, result.dir = NULL,
                           bedTools.dir, prefix = NULL, report = TRUE, prop.threshold = 0.001, delta = 0.1,
                           maf.control = 0.01, tau = 2, sigma = 0.1, len.threshold.K = 10,
                           group.length.threshold = 2, gain.threshold = log2(1.2), loss.threshold = log2(0.8),
                           WGD = 1.35, pos.prop.threhold = 0.6, pos.log2ratio.threhold = 0.75, Normalized = TRUE){

  if(is.null(result.dir)) result.dir <- "result"
  if(is.null(working.dir)) working.dir <- "working"

  len.threshold <- len.threshold.K

  library(inline)
  library(Rcpp)

  src <- '
    NumericVector X(X_);
    double tau = as<double>(tau_);
    NumericMatrix ans (X.size(), X.size());
    int ii, jj;
    tau = pow(tau, 2);
    for (ii=0; ii<X.size()-1; ii++){
      for (jj=ii+1; jj<X.size(); jj++){
        double sum = 0;
        sum = pow(X(ii) - X(jj), 2);
        double res = exp (- sum/tau);
        ans(ii, jj) = res;
        ans(jj, ii) = res;
      }
    }
    return(ans);'

  adj_vec_inline_Rcpp <- cxxfunction(signature(X_="numeric", tau_="numeric"), body = src, plugin="Rcpp")

  find.cluster <- function(signal.segments, tau = tau, sigma = sigma, k.range = c(1:7)){
    norm.likelihood <- function(x, grouping, sigma){
      group.mean <- tapply(x, grouping, mean)
      res <- sum(dnorm(x, group.mean[grouping], sd = sigma, log = T))
      return(res)
    }
    adj.segments <- adj_vec_inline_Rcpp(signal.segments, tau)
    dist.segments <- dist(adj.segments)
    hl.segments <- hclust(dist.segments, "average")
    likelihood <- NULL
    for(k in k.range){
      group.k <- try(cutree(hl.segments, k = k))
      if(class(group.k) != "try-error"){
        likelihood <- c(likelihood, norm.likelihood(signal.segments, group.k, sigma = sigma))
      } else {
        likelihood <- c(likelihood, -10)
      }
    }
    BIC <- -2*likelihood + c(k.range)*(log(length(signal.segments)) - log(2*pi))
    select.k <- k.range[which.min(BIC)]
    select.group <- cutree(hl.segments, k = select.k)
    return(select.group)
  }

  options(scipen = 50)

  if(class(Segment) != "Segment")
    stop("Invalid input class for purityEstimate().")

  segmentMethod <- Segment$segmentMethod

 if(!is.null(Segment$segmentNormalized) & Normalized == TRUE){
   segRes <- Segment$segmentNormalized
   Segment$WholeGenomeDoubling <- CallWholeGenomeDoubling(segRes, WGD = WGD, pos.prop.threhold = pos.prop.threhold, pos.log2ratio.threhold = pos.log2ratio.threhold)$WholeGenomeDoubling
   if(Segment$WholeGenomeDoubling == TRUE){
     segRes <- Segment$segmentUnadjusted
   }
 } else {
   segRes <- Segment$segmentUnadjusted
 }

  segRes$ratio <- 2^segRes[, "log2ratio"]
  subsegRes <- segRes[segRes[, "chr"] %in% c(1:22), ]
  segLens <- subsegRes[, "end"] - subsegRes[, "start"]
  ploidy <- sum(segLens*subsegRes[, "ratio"])/sum(segLens)

  subsegRes <- subsegRes[, c("chr", "start", "end", "ratio")]
  subsegRes <- subsegRes[segLens >= len.threshold & subsegRes[, "ratio"] < 2, ]
  a.sample.cluster <- find.cluster(subsegRes[, 4], tau = tau, sigma = sigma)
  num.clusters <- length(unique(a.sample.cluster))

  if(num.clusters > 1) {

    group.mean <- tapply(subsegRes[, 4], a.sample.cluster, mean)
    group.length <- tapply(subsegRes[, 4], a.sample.cluster, length)
    clustering <- data.frame(subsegRes, a.sample.cluster)
    write.table(clustering, paste0(working.dir, "/clustering.bed"), col.names = F, row.names = F, sep = "\t", quote = F)

    ff <- paste0(bedTools.dir, " -a ", working.dir, "/tumor.MAF.highcut.bed", " -b ", working.dir, "/clustering.bed -wa -wb > ",
               working.dir, "/tumor.MAF_clustering.bed")
    system(ff)

    maf.clustering <- read.delim(paste0(working.dir, "/tumor.MAF_clustering.bed"), header = F)
    maf <- tapply(maf.clustering[, 4], maf.clustering[, 9], median)
    if(length(maf) < length(group.mean)){
      group.mean <- group.mean[names(maf)]
      group.length <- group.length[names(maf)]
    }

    total.length <- tapply(c(maf.clustering[, 7] - maf.clustering[, 6]+1), maf.clustering[, 9], sum)
    total.prop <- total.length/sum(total.length)
    group.mean <- group.mean[total.prop > prop.threshold]
    group.length <- group.length[total.prop > prop.threshold]
    maf <-  maf[total.prop > prop.threshold]

    if(all(abs(diff(sort(maf))) <= maf.control)){
      admix <- NA
    } else {
      copy.neutral <- which.min(abs(group.mean-1))
      dis <- (group.mean[copy.neutral] - 1)
      group.mean <- group.mean - dis
      subsegRes[, 4] <- subsegRes[, 4] - dis
      group.label <- which(group.mean < 1-delta & group.length >= group.length.threshold)

      #### admixture: tumor
      tt <- copy.neutral == group.label
      if(length(tt) == 1) {
        tt_flag <- tt
      } else {
        tt_flag <- FALSE
      }
      if(length(group.label) == 0 | tt_flag) {
        group.label2 <- which(group.mean > 1+delta & group.length >= group.length.threshold)
        if(length(group.label2) > 0 & any(group.mean[group.label2] <  1.5)) {
          admix <- 2*min(group.mean[group.label2]) - 2
          if(admix > 1) admix <- 1
        } else { admix <- NA }
      } else if(length(group.label) > 1 & all(abs(diff(group.mean[group.label])) < delta)) {
        collapsed <- mean(subsegRes[a.sample.cluster%in%group.label, 4])
        admix <- 2 - 2*collapsed
      } else if(length(group.label) == 1) {
        admix <- 2 - 2*group.mean[group.label]
        if(admix > 1) admix <- 1
      } else if(length(group.label) > 1) {
        group.label <- group.label[!tt]
        candidate <- group.mean[group.label]
        candidate.onecopy <- candidate[candidate > 0.5]
        candidate.zerocopy <- candidate[candidate < 0.5]
        if(length(candidate.onecopy)>0){
          admix <- 2 - 2*candidate.onecopy[which.min(candidate.onecopy-0.5)]
        } else if(length(candidate.zerocopy)>0){
          admix <- 1 - min(candidate.zerocopy)/2
        }
      } else { admix <- NA }
    }

    purity <- admix

    if(!is.na(purity)){

      segRes <- Segment$segmentUnadjusted <- correctBypurity(Segment$segmentUnadjusted, purity)
      if(!is.null(Segment$segmentNormalized)){
        segRes <- Segment$segmentNormalized <- correctBypurity(Segment$segmentNormalized, purity)
      }

      res <- list(purity, ploidy*2)
      names(res) <- c("Purity", "Ploidy")
      Segment$PurityPloidy <- res

    } else {
      purity <- NA
      if(abs(ploidy*2-2) <= 0.05){

        res <- list(purity, ploidy*2)
        Segment$segmentUnadjusted <- correctBypurity(Segment$segmentUnadjusted, 1)
        if(!is.null(Segment$segmentNormalized)){
          Segment$segmentNormalized <- correctBypurity(Segment$segmentNormalized, 1)
        }

      } else {

        res <- list("<10%", ploidy*2)
        Segment$segmentUnadjusted <- correctBypurity(Segment$segmentUnadjusted, 0.1)
        if(!is.null(Segment$segmentNormalized)){
          Segment$segmentNormalized <- correctBypurity(Segment$segmentNormalized, 0.1)
        }

      }

      names(res) <- c("Purity", "Ploidy")
      Segment$PurityPloidy <- res
    }
  } else {
    purity <- NA
    if(abs(ploidy*2-2) <= 0.05){
      res <- list(purity, ploidy*2)
      Segment$segmentUnadjusted <- correctBypurity(Segment$segmentUnadjusted, 1)
      if(!is.null(Segment$segmentNormalized)){
        Segment$segmentNormalized <- correctBypurity(Segment$segmentNormalized, 1)
      }

    } else {
      res <- list("<10%", ploidy*2)

      Segment$segmentUnadjusted <- correctBypurity(Segment$segmentUnadjusted, 0.1)
      if(!is.null(Segment$segmentNormalized)){
        Segment$segmentNormalized <- correctBypurity(Segment$segmentNormalized, 0.1)
      }
    }
    names(res) <- c("Purity", "Ploidy")
    Segment$PurityPloidy <- res
  }

  segRes <- Segment$segmentUnadjusted
  subsegRes <- segRes[segRes[, "chr"] %in% c(1:22), ]
  segLens <- subsegRes[, "end"] - subsegRes[, "start"]
  ploidy <- sum(segLens*subsegRes[, "correctedRatiobyPurity"])/sum(segLens)
  Segment$PloidyUnadjusted <- ploidy*2

  if(!is.null(Segment$segmentNormalized)){
    segRes <- Segment$segmentNormalized
    subsegRes <- segRes[segRes[, "chr"] %in% c(1:22), ]
    segLens <- subsegRes[, "end"] - subsegRes[, "start"]
    ploidy <- sum(segLens*subsegRes[, "correctedRatiobyPurity"])/sum(segLens)
    Segment$PurityPloidy["Ploidy"] <- ploidy*2

  } else {
    Segment$PurityPloidy["Ploidy"] <- Segment$PloidyUnadjusted
  }

  if(Segment$PurityPloidy["Ploidy"] >= 3 | Segment$PloidyUnadjusted >= 2.7){
    Segment$WholeGenomeDoubling <- TRUE
  } else {
    Segment$WholeGenomeDoubling <- FALSE
  }

  if(!is.na(Segment$PurityPloidy["Purity"])){
    if(!is.na(purity)){
      if(purity < 0 ){
        purity <- NA
        res <- list(purity, ploidy*2)
        names(res) <- c("Purity", "Ploidy")
        Segment$PurityPloidy <- res
      }
    }
  }

#  if(Segment$PurityPloidy["Feature"] == "Normal Like") {
#    normallike <- TRUE
#  } else {
#    normallike <- FALSE
#  }

  reportname <- "SynthEx_sample_stats.bed"
  sampleinfo <- data.frame(Segment$PurityPloidy[c("Purity", "Ploidy")],
                           Segment$WholeGenomeDoubling, Segment$segmentMethod)

  colnames(sampleinfo) <- c("purity", "ploidy", "WholeGenomeDoubling", "segmentMethod")
  if(!is.null(prefix)){
    write.table(sampleinfo, paste0(result.dir, "/", prefix, "_", reportname), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  } else {
    write.table(sampleinfo, paste0(result.dir, "/", reportname), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }

  allsegRes <- Segment$segmentUnadjusted

  if(ncol(Segment$segmentUnadjusted) == 7){
    Event <- ifelse(allsegRes[, "correctedlog2ratiobyPurity"] >= gain.threshold, "Gain",
               ifelse(allsegRes[, "correctedlog2ratiobyPurity"] <= loss.threshold, "Loss", ""))
      Copy <- round(Segment$segmentNormalized[, "correctedRatiobyPurity"]*2, 2)
      IntegerCopy <- round(Copy)
      eventinfo <- data.frame(allsegRes[, c("chr", "start", "end", "log2ratio", "correctedlog2ratiobyPurity")], Event)
      copynumber <- data.frame(Segment$segmentNormalized[, c("chr", "start", "end", "log2ratio", "correctedlog2ratiobyPurity")], Copy, IntegerCopy)
      colnames(eventinfo) <- c("chr", "start", "end", "log2ratio", "correctedlog2ratiobyPurity", "event")
      colnames(copynumber) <- c("chr", "start", "end", "normalizedlog2ratio", "normalizedCorrectedlog2ratiobyPurity", "Copy", "integerCopy")
  } else {
      Event <- ifelse(allsegRes[, "log2ratio"] >= gain.threshold, "Gain",
                        ifelse(allsegRes[, "log2ratio"] <= loss.threshold, "Loss", ""))
      Copy <- round(Segment$segmentNormalized[, "ratio"]*2, 2)
      IntegerCopy <- round(Copy)
      eventinfo <- data.frame(allsegRes[, c("chr", "start", "end", "log2ratio")], Event)
      copynumber <- data.frame(Segment$segmentNormalized[, c("chr", "start", "end", "log2ratio")], Copy, IntegerCopy)
      colnames(eventinfo) <- c("chr", "start", "end", "log2ratio", "event")
      colnames(copynumber) <- c("chr", "start", "end", "normalizedlog2ratio", "Copy", "integerCopy")
    }

  if(!is.null(prefix)){
    write.table(eventinfo, paste0(result.dir, "/", prefix, "_Event.bed"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    write.table(copynumber, paste0(result.dir, "/", prefix, "_Copynumber.bed"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  } else {
    write.table(eventinfo, paste0(result.dir, "/Event.bed"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    write.table(copynumber, paste0(result.dir, "/Copynumber.bed"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  Segment$sampleinfo <- sampleinfo
  Segment$event <- eventinfo
  Segment$copynumber <- copynumber
  return(Segment)

}



