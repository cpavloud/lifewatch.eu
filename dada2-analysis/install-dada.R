#!/usr/bin/env Rscript

source('http://bioconductor.org/biocLite.R')
biocLite('ShortRead', suppressUpdates = TRUE)
biocLite('dada2', suppressUpdates = TRUE)
install.packages('ggplot2', repos='https://cloud.r-project.org/')
install.packages('igraph', repos='https://cloud.r-project.org/')
install.packages('reshape', repos='https://cloud.r-project.org/')
biocLite('phyloseq', suppressUpdates = TRUE)
install.packages('argparser')
#install.packages('optparse')
