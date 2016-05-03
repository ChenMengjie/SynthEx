SynthExUserGuide <- function(view=TRUE)
  #	Written by Mengjie, Apr 2016; modified from edgeRUsersGuide() in package edgeR.
  #	Find User's Guide
{
  f <- system.file("doc", "SynthExUserGuide.pdf", package="SynthEx")
  if(view) {
    if(.Platform$OS.type == "windows")
      shell.exec(f)
    else
      system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
  }
  return(f)
}
