rm(list=ls())

library(dplyr)
library(glmnet)

panel = read.table("panel_chrom1.txt", header=F)
wsaf = read.table("PG0402-C_chrom1.wsaf", header=F)$V1
myfit = glmnet( x = as.matrix(panel),
                y = wsaf,
                lambda = 1/seq(2, 100, length.out = 100),
                type.gaussian = "naive",
                lower.limits= 0)


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


######################### Normalizing y #######################
y = wsaf
w = rep(1, length(wsaf))
w=w/sum(w)
v=sqrt(w)

ym=sum(w*y)
y=v*(y-ym)
ys=sqrt(sum(y*y))
y=y/ys

####### my calculation
ym = mean(wsaf)
ytmp = (wsaf - ym)
ytmp / sqrt(sum(ytmp^2))


########################## Normalizing x #####################


#10710 do 10711 j=1,ni
#            if(ju(j).eq.0)goto 10711
#            xm(j)=dot_product(w,x(:,j))
#            x(:,j)=v*(x(:,j)-xm(j))
#            xv(j)=dot_product(x(:,j),x(:,j))
#            if(isd.gt.0) xs(j)=sqrt(xv(j))
#      10711 continue
#      10712 continue

#      if(isd .ne. 0)goto 10731
#            xs=1.0
#            goto 10741
#      10731 continue

#      10750 do 10751 j=1,ni
#            if(ju(j).eq.0)goto 10751
#            x(:,j)=x(:,j)/xs(j)
#      10751 continue
#      10752 continue
#            xv=1.0
#      10741 continue
#10721 continue
