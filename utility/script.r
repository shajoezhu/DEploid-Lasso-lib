#!/usr/bin/env Rscript
rm(list=ls())
library(dplyr)

args = (commandArgs(TRUE))

inFileName = args[1]

system(paste("grep TABLE ", inFileName, " > ", inFileName, ".TABLE", sep=""))
system(paste("grep BETA ", inFileName, " > ", inFileName, ".BETA", sep=""))

mytable = read.table(paste(inFileName, ".TABLE", sep=""), header=T, stringsAsFactors = F)[,-1]
mybeta = read.table(paste(inFileName, ".BETA", sep=""), header=T, stringsAsFactors = F)[,-1]


pdf("tmp.pdf")

plot(mytable$df)

plot(mytable$rsq)

plot(diff(mytable$rsq))

plot(c(0, max(mytable$L1norm)),
     c(min(mybeta[,-1]), max(mybeta[,-1])), type="n")
for (i in 1:dim(mybeta)[1]){
    lines(c(0,mytable$L1norm), c(0,mybeta[i,-1]))
}

dev.off()
