fun.divide.to.seg <- function(hapLength, by = 50){
#    return(floor(seq(1, hapLength, length.out =numSeg)))
    myseq = seq(1, hapLength, by = by)
    if ((hapLength-1)%%by != 0){
        myseq = c(myseq, hapLength)
    }
    return(myseq)
}

getIndex3 <- function (mixedSample, ref1, ref2, ref3, prop){
#1	1	2	2	3	3
#2	3	1	3	2	1
#3	2	3	1	1	2
#    prop = c(.25,.25,.5)
    prop = c(1,1,1)
    mySums = c ( sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
                 sum(mixedSample[,1] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] ,
                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,3] != ref3)*prop[3] ,
                 sum(mixedSample[,2] != ref1)*prop[1] + sum(mixedSample[,3] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,2] != ref2)*prop[2] + sum(mixedSample[,1] != ref3)*prop[3] ,
                 sum(mixedSample[,3] != ref1)*prop[1] + sum(mixedSample[,1] != ref2)*prop[2] + sum(mixedSample[,2] != ref3)*prop[3] )
#print(mySums)
#print("")
    case = which.min(mySums)
    if ( case == 1 ){
        return (c(1,2,3))
    } else if ( case == 2 ){
        return (c(1,3,2))
    } else if ( case == 3 ){
        return (c(2,1,3))
    } else if ( case == 4 ){
        return (c(2,3,1))
    } else if ( case == 5 ){
        return (c(3,2,1))
    } else if ( case == 6 ){
        return (c(3,1,2))
    }
}


getIndex2 <- function (mixedSample, ref1, ref2){
    sum1 = sum(mixedSample[,1] != ref1) + sum(mixedSample[,2] != ref2)
    sum2 = sum(mixedSample[,2] != ref1) + sum(mixedSample[,1] != ref2)
    if ( sum1 < sum2 ){
        return (c(1,2))
    } else {
        return (c(2,1))
    }
}


exportSites <- function ( CHROM, POS, tmpIndex, logFileName){
    for ( index in tmpIndex ){
        cat ( as.character(CHROM[index]), ",", POS[index], "\n", file = logFileName, append = T)
    }
}


fun.computeErrors2 <- function(hap, ref1, ref2){
    switchError = c(0, 0)
    mutError = c(0, 0)

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    if ( nhap == 2 ){

        index.of.seg = fun.divide.to.seg(haplength)
        countSwitch = matrix(NA, ncol = length(index.of.seg)-1, nrow = 2)

        for ( i in 1:(length(index.of.seg)-1) ){
            tmpIndex = c(index.of.seg[i]:(index.of.seg[i+1]-1))
            # k = 1
            if ( sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref2[tmpIndex]) ){
                hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
                countSwitch[1,i] = 1
            } else {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref2[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 7
                countSwitch[1,i] = 2
            }
#            } else if ( sum(hap[tmpIndex,1] != ref3[tmpIndex]) <= sum(hap[tmpIndex,1] != ref2[tmpIndex]) ) {
#                hap[tmpIndex,1][hap[tmpIndex,1] != ref3[tmpIndex]] = 2
##                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 8
#                countSwitch[1,i] = 3
#            }

            # k = 2
            if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) ){
                hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                countSwitch[2,i] = 2
            } else {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
#                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 8
                countSwitch[2,i] = 1
            }
#            } else if ( sum(hap[tmpIndex,2] != ref1[tmpIndex]) <= sum(hap[tmpIndex,2] != ref3[tmpIndex]) ) {
#                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
#                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
#                countSwitch[2,i] = 1
#            }
        }
        strain = cbind(c(NA,NA), countSwitch)
        strainNext = cbind(countSwitch, c(NA,NA))
        switchError = rowSums(strain != strainNext, na.rm=T)
        mutError = colSums(hap==2)
    }

    return ( list ( hap = hap,
                    switchError = switchError,
                    mutError = mutError ))

}


fun.plotHapWithProp <- function( hap, prop, fig.title, max.at ){

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    xrange = c(0, max.at)
    yrange = c(0, 1)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 4, cex.axis=2)
#    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 3, cex.axis=2, xaxt="n")

    plot( xrange, yrange, type= "n", xlim=xrange, ylim = yrange, ylab="", main=fig.title, xlab = "", cex.lab = 2, cex.main = 1.4, cex.axis=2, xaxt="n")
    xleft = 0:(haplength-1)
    xright = xleft+1
#    ycum = as.numeric(c(0, prop[1], 1))
    ycum = as.numeric(c(0, cumsum(as.numeric(prop))))
#cat(ycum, "\n")

    for ( k in c(1:nhap) ){
      tmpHap = hap[,k]
      ybottom = ycum[k]
      ytop = ycum[k+1]
      rect(xleft, ybottom, xright, ytop, col = tmpHap , border = "transparent")
    }
    lab.at = c(1, 400, 800, 1200, 1600, 2000, 2369)
    axis(1, at=lab.at, labels=lab.at, las=1, lwd = 0, cex=2, cex.axis=2.4)

}

fun.computeErrors3 <- function(hap, ref1, ref2, ref3){
    switchError = c(0, 0, 0)
    mutError = c(0, 0, 0)

    haplength = dim(hap)[1]
    nhap = dim(hap)[2]

    if ( nhap == 3 ){
#        hap[hap[,1] != ref1, 1 ] = 2
#        hap[hap[,2] != ref2, 2 ] = 2
#        hap[hap[,3] != ref3, 3 ] = 2

        index.of.seg = fun.divide.to.seg(haplength)
        countSwitch = matrix(NA, ncol = length(index.of.seg)-1, nrow = 3)

        for ( i in 1:(length(index.of.seg)-1) ){
            tmpIndex = c(index.of.seg[i]:(index.of.seg[i+1]-1))
            # k = 1
            if ( sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref2[tmpIndex]) & sum(hap[tmpIndex,1] != ref1[tmpIndex]) < sum(hap[tmpIndex,1] != ref3[tmpIndex]) ){
                hap[tmpIndex,1][hap[tmpIndex,1] != ref1[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 5
                hap[tmpIndex,1][hap[tmpIndex,1] == 1] = 5 # turn off the black color
                countSwitch[1,i] = 1
            } else if ( sum(hap[tmpIndex,1] != ref2[tmpIndex]) < sum(hap[tmpIndex,1] != ref3[tmpIndex]) ) {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref2[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 0] = 7
                hap[tmpIndex,1][hap[tmpIndex,1] == 1] = 7 # turn off the black color
                countSwitch[1,i] = 2
            } else if ( sum(hap[tmpIndex,1] != ref3[tmpIndex]) <= sum(hap[tmpIndex,1] != ref2[tmpIndex]) ) {
                hap[tmpIndex,1][hap[tmpIndex,1] != ref3[tmpIndex]] = 2
                hap[tmpIndex,1][hap[tmpIndex,1] == 1] = 0 # turn off the black color
                countSwitch[1,i] = 3
            }

            # k = 2
            if ( sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) & sum(hap[tmpIndex,2] != ref2[tmpIndex]) < sum(hap[tmpIndex,2] != ref3[tmpIndex]) ){
                hap[tmpIndex,2][hap[tmpIndex,2] != ref2[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 7
                hap[tmpIndex,2][hap[tmpIndex,2] == 1] = 7 # turn off the black color
                countSwitch[2,i] = 2
            } else if ( sum(hap[tmpIndex,2] != ref3[tmpIndex]) < sum(hap[tmpIndex,2] != ref1[tmpIndex]) ) {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref3[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 1] = 0 # turn off the black color
                countSwitch[2,i] = 3
            } else if ( sum(hap[tmpIndex,2] != ref1[tmpIndex]) <= sum(hap[tmpIndex,2] != ref3[tmpIndex]) ) {
                hap[tmpIndex,2][hap[tmpIndex,2] != ref1[tmpIndex]] = 2
                hap[tmpIndex,2][hap[tmpIndex,2] == 0] = 5
                hap[tmpIndex,2][hap[tmpIndex,2] == 1] = 5 # turn off the black color
                countSwitch[2,i] = 1
            }

            # k = 3
            if ( sum(hap[tmpIndex,3] != ref3[tmpIndex]) < sum(hap[tmpIndex,3] != ref2[tmpIndex]) & sum(hap[tmpIndex,3] != ref3[tmpIndex]) < sum(hap[tmpIndex,3] != ref1[tmpIndex]) ){
                hap[tmpIndex,3][hap[tmpIndex,3] != ref3[tmpIndex]] = 2
                hap[tmpIndex,3][hap[tmpIndex,3] == 1] = 0 # turn off the black color
                countSwitch[3,i] = 3
            } else if ( sum(hap[tmpIndex,3] != ref2[tmpIndex]) < sum(hap[tmpIndex,3] != ref1[tmpIndex]) ) {
                hap[tmpIndex,3][hap[tmpIndex,3] != ref2[tmpIndex]] = 2
                hap[tmpIndex,3][hap[tmpIndex,3] == 0] = 7
                hap[tmpIndex,3][hap[tmpIndex,3] == 1] = 7 # turn off the black color
                countSwitch[3,i] = 2
            } else if ( sum(hap[tmpIndex,3] != ref1[tmpIndex]) <= sum(hap[tmpIndex,3] != ref2[tmpIndex]) ) {
                hap[tmpIndex,3][hap[tmpIndex,3] != ref1[tmpIndex]] = 2
                hap[tmpIndex,3][hap[tmpIndex,3] == 0] = 5
                hap[tmpIndex,3][hap[tmpIndex,3] == 1] = 5 # turn off the black color
                countSwitch[3,i] = 1
            }
        }
        strain = cbind(c(NA,NA,NA), countSwitch)
        strainNext = cbind(countSwitch, c(NA,NA,NA))
        switchError = rowSums(strain != strainNext, na.rm=T)
        mutError = colSums(hap==2)
    }

    return ( list ( hap = hap,
                    switchError = switchError,
                    mutError = mutError ))

}


measure.error.joe<-function(h.pair, h.pair.true, rel.cost.switch=2, do.plot=FALSE) {
    l <- ncol(h.pair);
    n.hap <- nrow(h.pair)
    possible.permn = combinat::permn(1:n.hap)
    n.permn = length(possible.permn)
    v<-rep(0, n.permn);
    vn<-v;

    tb<-array(0, c(n.permn, l));

    for ( j in 1:n.permn){
        v[j] = sum(h.pair[,1]!=h.pair.true[possible.permn[[j]],1]);
    }

    ee <- rep(0, n.permn)
    for (i in 2:l) {
        for ( j in 1:n.permn){
            ee[j] = sum(h.pair[,i]!=h.pair.true[possible.permn[[j]],i]);
            ones = rep(1, n.permn)
            ones[j] = 0
            tmp <- v + rel.cost.switch * ones
            vn[j] <- min(tmp) + ee[j]
            tb[j, i] <- which.min(tmp)
        }
        v<-vn;
    }

    #decode
    wm<-which.min(v);
    op<-array(0, c(1,l));
    n.s<-0;

    if (wm!=0){
        n.gt<-sum(h.pair[,l]!=h.pair.true[possible.permn[[wm]],l]);
    }


    op[l]<-wm;
    for (i in l:2) {
        wmp<-tb[wm,i];
        if (wmp!=wm) n.s<-n.s+1;

        if (wmp!=0){
            n.gt <- n.gt + sum(h.pair[,i-1] != h.pair.true[possible.permn[[wmp]],i-1]);
        }

        wm<-wmp;
        op[i-1]<-wm;
    }
	cat("\nDecoding gives:\nNo. switches:\t", n.s, "\nNo. GT errs:\t", n.gt, "\n\n");

	if (do.plot) {
        par(mfrow=c(2,1))
		plot(0,0,type="n", xlab="Position", ylab="", yaxt="n", xlim=c(0,ncol(h.pair)), bty="n", ylim=c(-0.5,1.5));
		image(x=1:ncol(h.pair), z=t(rbind(h.pair, rep(0, ncol(h.pair)), h.pair.true)), col=c("white", "black"),
			add=T);
		del<-which(diff(op[1,])!=0);
		if (length(del)>0) points(x=del, y=rep(0.5, length(del)), pch="|", col="red");

        for ( j in 1:n.hap ){
            wd1<-which(h.pair[j,] != h.pair.true[j,]);
    #
            if (length(wd1)>0) {
                points(x=wd1, y=rep(1.2, length(wd1)), pch=25, col=j+1, cex=0.5);
            }

        }
        image(t(op))
    }

#	return(c(n.s, n.gt));
        return (list(switchError = n.s,
                 mutError = n.gt,
                 op = op) )
}

measure.error.joe.2<-function(h.pair, h.pair.true, rel.cost.switch=2, do.plot=FALSE) {
    hapAndError = measure.error.joe(h.pair, h.pair.true, rel.cost.switch, do.plot)

    n.hap = nrow(h.pair)
#    cat("n.hap = ", n.hap, "\n")
    switchError = rep(0, n.hap)
    mutError = rep(0, n.hap)
    possible.permn = combinat::permn(1:n.hap)

    l = ncol(h.pair.true)
    hap = array(0, c(nrow(h.pair), l));
    for ( i in 1:l ){
        hap[,i] = possible.permn[[hapAndError$op[i]]]
    }

    for ( j in 1:n.hap ){
        switchError[j] = sum((hap[j,-l] - hap[j, -1]) != 0)
    }

    for ( i in 1:l ){
        hap[h.pair[,i] != h.pair.true[hap[,i],i], i ] = 0
    }

    for ( j in 1:n.hap ){
        mutError[j] = sum(hap[j,] == 0)
    }

    return ( list ( hap = t(hap),
                    switchError = switchError,
                    mutError = mutError ))

}

measure.error.joe.with.drop<-function(h.pair, h.pair.true, rel.cost.switch=2, rel.cost.drop = 1, do.plot=FALSE) {
    l <- ncol(h.pair);
    n.hap <- nrow(h.pair)
    possible.permn = combinat::permn(1:n.hap)
    current.len = length(possible.permn)

possible.permn[[7]] = c(1,2,1)
possible.permn[[8]] = c(1,1,2)
possible.permn[[9]] = c(2,1,1)

possible.permn[[10]] = c(1,3,1)
possible.permn[[11]] = c(1,1,3)
possible.permn[[12]] = c(3,1,1)

possible.permn[[13]] = c(2,3,2)
possible.permn[[14]] = c(2,2,3)
possible.permn[[15]] = c(3,2,2)

possible.permn[[16]] = c(2,1,2)
possible.permn[[17]] = c(2,2,1)
possible.permn[[18]] = c(1,2,2)

possible.permn[[19]] = c(3,1,3)
possible.permn[[20]] = c(3,3,1)
possible.permn[[21]] = c(1,3,3)

possible.permn[[22]] = c(3,2,3)
possible.permn[[23]] = c(3,3,2)
possible.permn[[24]] = c(2,3,3)

    n.permn = length(possible.permn)
print(n.permn)
    v<-rep(0, n.permn);
    vn<-v;

    tb<-array(0, c(n.permn, l));

    for ( j in 1:n.permn){
        v[j] = sum(h.pair[,1]!=h.pair.true[possible.permn[[j]],1]);
    }

    ee <- rep(0, n.permn)
    for (i in 2:l) {
        for ( j in 1:n.permn){
            ee[j] = sum(h.pair[,i]!=h.pair.true[possible.permn[[j]],i]);
            ones = rep(1, n.permn)
            ones[j] = 0
            drop.ones = rep(1, n.permn)
            drop.ones[j] = 0
            if (j<7){
                drop.ones[1:6] = 0
            } else {
                drop.ones[7:24] = 0
            }
            tmp <- v + rel.cost.switch * ones + rel.cost.drop*drop.ones
            vn[j] <- min(tmp) + ee[j]
            tb[j, i] <- which.min(tmp)
        }
        v<-vn;
    }

    #decode
    wm<-which.min(v);
    op<-array(0, c(1,l));
    n.s<-0;

    if (wm!=0){
        n.gt<-sum(h.pair[,l]!=h.pair.true[possible.permn[[wm]],l]);
    }


    op[l]<-wm;
    for (i in l:2) {
        wmp<-tb[wm,i];
        if (wmp!=wm) n.s<-n.s+1;

        if (wmp!=0){
            n.gt <- n.gt + sum(h.pair[,i-1] != h.pair.true[possible.permn[[wmp]],i-1]);
        }

        wm<-wmp;
        op[i-1]<-wm;
    }

	if (do.plot) {
        par(mfrow=c(2,1))
		plot(0,0,type="n", xlab="Position", ylab="", yaxt="n", xlim=c(0,ncol(h.pair)), bty="n", ylim=c(-0.5,1.5));
		image(x=1:ncol(h.pair), z=t(rbind(h.pair, rep(0, ncol(h.pair)), h.pair.true)), col=c("white", "black"),
			add=T);
		del<-which(diff(op[1,])!=0);
		if (length(del)>0) points(x=del, y=rep(0.5, length(del)), pch="|", col="red");

        for ( j in 1:n.hap ){
            wd1<-which(h.pair[j,] != h.pair.true[j,]);
    #
            if (length(wd1)>0) {
                points(x=wd1, y=rep(1.2, length(wd1)), pch=25, col=j+1, cex=0.5);
            }

        }
        image(t(op))
    }

    k.eff.permn = rep(0, n.permn)
    for (j in 1:n.permn){
        k.eff.permn[j] = length(unique(possible.permn[[j]]))
    }
    k.eff = k.eff.permn[op]

    drop.strain.permn = rep(0, n.permn)
    for (j in 1:n.permn){
        dropped = which(!c(1,2,3) %in% possible.permn[[j]])
        if ( length(dropped) > 0 ){
            drop.strain.permn[j] = dropped
        }
    }
    drop.strain =  drop.strain.permn[op]

    dropTimes = sum(diff(drop.strain) != 0)
    dropError = sum(drop.strain != 0)
    cat("\nDecoding gives:\nNo. switches:\t", n.s-dropTimes, "\nNo. GT errs:\t", n.gt, "\nNo. Drop errs:\t", dropError,"\n");
        return (list(switchError = n.s - dropTimes,
                 mutError = n.gt,
                 dropError = dropError,
                 op = op, drop.strain = drop.strain,
                 possible.permn = possible.permn) )
}


measure.error.joe.with.drop.2strain <-function(h.pair, h.pair.true, rel.cost.switch=2, do.plot=FALSE) {
    l <- ncol(h.pair);
    n.hap <- nrow(h.pair)
    possible.permn = combinat::permn(1:n.hap)

possible.permn[[3]] = c(1,1)
possible.permn[[4]] = c(2,2)

#possible.permn = list()
#possible.permn[[1]] = c(1,1)
#possible.permn[[2]] = c(2,2)
#possible.permn[[3]] = c(1,2)
#possible.permn[[4]] = c(2,1)

    n.permn = length(possible.permn)
    v<-rep(0, n.permn);
    vn<-v;

    tb<-array(0, c(n.permn, l));

    for ( j in 1:n.permn){
        v[j] = sum(h.pair[,1]!=h.pair.true[possible.permn[[j]],1]);
    }

#rel.cost.drop.current = 0
    ee <- rep(0, n.permn)
    same.path <- rep(0, n.permn)

    for (i in 2:l) {
#        rel.cost.drop.current = rel.cost.drop.current + 1
        for ( j in 1:n.permn){
            ee[j] = sum(h.pair[,i]!=h.pair.true[possible.permn[[j]],i]);

            ones = rep(1, n.permn)
            ones[j] = 0
            drop.ones = rep(1, n.permn)
            drop.ones[j] = 0
            if (j<3){
                drop.ones[1:2] = 0
            } else {
                drop.ones[3:4] = 0
            }
#            tmp <- v + rel.cost.switch * ones
            tmp <- v + rel.cost.switch * ones + same.path *drop.ones

            vn[j] <- min(tmp) + ee[j]
            transit.to = which.min(tmp)

            if ( j == transit.to ){
                same.path = same.path + 1
            } else {
                same.path = 0
            }

            tb[j, i] <- transit.to
        }


#        if ( rel.cost.drop.current > rel.cost.drop){
#            rel.cost.drop.current = 0
#        }
        v<-vn;
#        cat ("site ", i, " cost: ", v,"\n")
    }

    #decode
    wm<-which.min(v);
    op<-array(0, c(1,l));
    n.s<-0;

    if (wm!=0){
        n.gt<-sum(h.pair[,l]!=h.pair.true[possible.permn[[wm]],l]);
    }


    op[l]<-wm;
    for (i in l:2) {
        wmp<-tb[wm,i];
        if (wmp!=wm) n.s<-n.s+1;

        if (wmp!=0){
            n.gt <- n.gt + sum(h.pair[,i-1] != h.pair.true[possible.permn[[wmp]],i-1]);
        }

        wm<-wmp;
        op[i-1]<-wm;
    }

    k.eff.permn = rep(0, n.permn)
    for (j in 1:n.permn){
        k.eff.permn[j] = length(unique(possible.permn[[j]]))
    }
    k.eff = k.eff.permn[op]

    changed.permn.at = which(diff(op[1,]) != 0)
    changed.k.eff.at = which(diff(k.eff) != 0)

#print("changed.permn.at")
#print(changed.permn.at)
#print("changed.k.eff.at")
#print(changed.k.eff.at)


    drop.strain.permn = rep(0, n.permn)
    for (j in 1:n.permn){
        dropped = which(!c(1,2) %in% possible.permn[[j]])
        if ( length(dropped) > 0 ){
            drop.strain.permn[j] = dropped
        }
    }
    drop.strain =  drop.strain.permn[op]

#    definite.switch = sum(!changed.permn.at %in% changed.k.eff.at)
    dropTimes = sum(diff(drop.strain) != 0) #length(changed.k.eff.at)
    dropError = sum(drop.strain != 0)
    cat("\nDecoding gives:\nNo. switches:\t", n.s - dropTimes, "\nNo. GT errs:\t", n.gt, "\nNo. Drop switches", dropTimes, "\nNo. Drop errs:\t", dropError,"\n");


    hap.pair.true.idx = array(unlist(possible.permn[op]), c(2, length(op)))
    genotype.error = rep(0, n.hap)
    switch.error = rep(0, n.hap)
    drop.error = rep(0, n.hap)
    for ( j in 1:n.hap ){
        wd1 = c()
        for ( i in 1:l ){
            if ( h.pair[j,i] != h.pair.true[hap.pair.true.idx[j,i],i] ){
                wd1 = c(wd1, i)
            }
        }
        genotype.error[j] = length(wd1)
    }

    switch.error[1] = (sum(diff(hap.pair.true.idx[1,])!=0))
    switch.error[2] = (sum(diff(hap.pair.true.idx[2,])!=0))
    drop.error[1] = (sum(drop.strain == 1))
    drop.error[2] = (sum(drop.strain == 2))

	if (do.plot) {
#        par(mfrow=c(2,1))
        layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = T))
		plot(0,0,type="n", xlab="Position", ylab="", yaxt="n", xlim=c(0,ncol(h.pair)), bty="n", ylim=c(-0.5,1.5));
		image(x=1:ncol(h.pair), z=t(rbind(h.pair, rep(0, ncol(h.pair)), h.pair.true)), col=c("white", "black"),
			add=T);
		del<-which(diff(op[1,])!=0);
		if (length(del)>0) {points(x=del, y=rep(0.5, length(del)), pch="|", col="red");}


        for ( j in 1:n.hap ){
            wd1 = c()
            for ( i in 1:l ){
                if ( h.pair[j,i] != h.pair.true[hap.pair.true.idx[j,i],i] ){
                    wd1 = c(wd1, i)
                }
            }

            if (length(wd1)>0) {
                points(x=wd1, y=rep(1.2, length(wd1)), pch=25, col=j+1, cex=0.5);
#                print(length(wd1))
            }

        }
		plot(0,0,type="n", xlab="Position", ylab="", yaxt="n", xlim=c(0,ncol(h.pair)), bty="n", ylim=c(-0.5,1.5));
#        image(x=1:ncol(h.pair), z = t(rbind(op, array(drop.strain,c(1,l)))), col = heat.colors(4), add=T)
        image(x=1:ncol(h.pair), z = t(op), col = heat.colors(4), add=T)
    }


        return (list(switchError = n.s - dropTimes,
                 mutError = n.gt,
                 dropError = dropError,
                 op = op, drop.strain = drop.strain,
                 possible.permn = possible.permn, tb = tb,
                 switch.error = switch.error,
                 drop.error = drop.error,
                 genotype.error = genotype.error )
                 )
}


measure.error.joe.with.drop.2strain.2<-function(h.pair, h.pair.true, rel.cost.switch=2, do.plot=FALSE) {
    hapAndError = measure.error.joe.with.drop.2strain(h.pair, h.pair.true, rel.cost.switch, do.plot)

    n.hap = nrow(h.pair)
#    cat("n.hap = ", n.hap, "\n")
    switchError = rep(0, n.hap)
    mutError = rep(0, n.hap)
    dropTimes = rep(0, n.hap)
    dropError = rep(0, n.hap)

#print(hapAndError$possible.permn)
    l = ncol(h.pair.true)
    hap = array(0, c(nrow(h.pair), l));
    for ( i in 1:l ){
        hap[,i] = hapAndError$possible.permn[[hapAndError$op[i]]]
    }
#        print(hapAndError$drop.strain)
     for ( j in 1:n.hap ){
       dropTimes[j] = sum(diff(hapAndError$drop.strain == j) != 0)
    }

#    for ( j in 1:n.hap ){
#        switchError[j] = sum(diff(hapAndError$op[1,]) != 0) - sum(dropTimes)
#        if (switchError[j] < 0) switchError[j] = 0
#    }

    for ( j in 1:n.hap ){
        switchError[j] = sum(diff(hap[j,]) != 0) - dropTimes[j]
        if (switchError[j] < 0) switchError[j] = 0
    }

    for ( j in 1:n.hap ){
        dropError[j] = sum(hapAndError$drop.strain == j)
    }

    for ( i in 1:l ){
        hap[h.pair[,i] != h.pair.true[hap[,i],i], i ] = 0
    }

    for ( j in 1:n.hap ){
        mutError[j] = sum(hap[j,] == 0)
    }

    cat("switchError ", switchError, "\n")
    cat("mutError ", mutError,"\n")
    cat("dropError ", dropError,"\n")

    return ( list ( hap = t(hap),
                    switchError = switchError,
                    mutError = mutError,
                    dropError = dropError ))

}
