Lookup <- function(x, show = 6){
  if(class(x) == "list"){
    len <- length(x)
    varnames <- names(x)
    for(i in 1:len){
      print(paste0("No.", i," element is ", varnames[i], ":"))
      if(class(x[[i]]) == "data.frame" | class(x[[i]]) == "matrix"){
        if(ncol(x[[i]]) > show){
          print(head(x[[i]][, 1:show]))
        } else {
          print(head(x[[i]]))
        }
      }
      if(class(x[[i]]) == "vector" | class(x[[i]]) == "character") print(head(x[[i]]))
      if(class(x[[i]]) == "list") print(summary(x[[i]]))
    }
  } else {
    print(summary(x))
  }
}
