rm(list=ls())

library(dplyr)
library(glmnet)

panel = read.table("panel_chrom1.txt", header=F)
wsaf = read.table("PG0402-C_chrom1.wsaf", header=F)$V1
myfit = glmnet(x=as.matrix(panel), y = wsaf, lambda = 1/seq(2, 100, length.out = 100))


lambdas = 1/seq(2, 100, length.out = 100)
results = list()
for ( i in 1:100 ){
    results[[i]] = admm.lasso(as.matrix(cbind(0,panel)),wsaf, lambdas[i])
}

a = admm.lasso(as.matrix(cbind(0,panel)),wsaf, 0.5)
# first index is for the intercept


which(results[100][[1]]$x > 0.01)
 [1] 13 21 25 27 30 38 62 69 76 88 92
> which(myfit$beta[,100] > 0.01)
V13 V21 V27 V62 V67 V88 V92
 13  21  27  62  67  88  92


y = wsaf
ybar = mean(y)
nulldev = sum((y - ybar)^2)

idx = 100
ytmp = myfit$a0[idx] + as.matrix(panel) %*% as.numeric(myfit$beta[,idx])

1 - sum((y - ytmp)^2)/nulldev

myfit$dev.ratio[idx]


ytmp = as.matrix(cbind(0,panel)) %*% a$x
1 - sum((y - ytmp)^2)/nulldev
