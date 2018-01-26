rm(list=ls())

library(dplyr)
library(glmnet)

panel = read.table("panel_chrom1.txt", header=F)
wsaf = read.table("PG0402-C_chrom1.wsaf", header=F)$V1
myfit = glmnet(x=as.matrix(panel), y = wsaf, lambda = 1/seq(2, 100, length.out = 100))
