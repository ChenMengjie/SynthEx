denoise <- function(data, k = 30, t = 2){
  #### Written by Mengjie Chen, Nov 2012
  #### smoothing data by replacing outliers (defined by deviation from mean in a sliding window with k points) with median  
  n <- 1:length(data)
  res <- tapply(n, 1:length(n), function(x){
    if(x<=k) return(data[c(1:(2*k))])
    else if(x>=length(n)-k) return(data[c((length(n)-2*k+1):length(n))])
    else return(data[c((x-k):(x+k-1))])
  })
  res <- matrix(unlist(res), byrow=T, ncol = 2*k)
  aa <- apply(res, 1, function(x){ mean(x, na.rm=T)})
  cc <- apply(res, 1, function(x){ median(x, na.rm=T)})
  bb <- apply(res, 1, function(x){ sd(x, na.rm=T)})
  flag <- (abs(data-aa) < t*bb) & !is.na(data)
  y <- data
  y[!flag] <- cc[!flag]
  return(y)
}
