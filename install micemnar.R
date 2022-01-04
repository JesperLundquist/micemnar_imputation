##Installerar micemnar från R:s arkiv##

url<-"https://cran.r-project.org/src/contrib/Archive/miceMNAR/miceMNAR_1.0.2.tar.gz"
pkgFile <- "miceMNAR_1.0.2.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(c("stats", "mvtnorm", "pbivnorm",  "GJRM",  "sampleSelection"))
install.packages(pkgs=pkgFile, type="source", repos=NULL)
unlink(pkgFile)
library(miceMNAR)
