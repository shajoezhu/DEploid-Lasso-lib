rm(list=ls())
library(lattice)

#pdf("switchVsMisCopyErrlogHi.pdf", width=8, height = 8)
pdf("HapErrors.withDrop.lasso.pdf", width=8, height = 12)
par(mar=c(4.1,5.1,2.1,2.1))

totalSites = 17530

simpleMix = paste("PG03", 89:94, "-C", sep="")
#mixedSampleNames = c(paste("PG03", 89:99, "-C", sep=""), paste("PG040", 0:9, "-C", sep=""), paste("PG04", 10:15, "-C", sep=""))

#plot(c(0,1), c(0, 2000), type="n")
mytable = data.frame ( prop = numeric(0), switchError = numeric(0), missCopy = numeric(0), dropError = numeric(0), sampleName = character(), strainName = character(), labStrain = character(), marker = character(), color = character(), stringsAsFactors = F)
mytablei = 1


#metaData = read.table("pf3k_release_5_metadata.txt", header = T, comment.char="", sep="\t", stringsAsFactors=F)

sample2 = read.table("~/dEploidPaper/validation/labSampleNames2Strains", header=F, stringsAsFactors=F)$V1
for (samplei in sample2){
    for ( seed in 1:15 ){
        prefix = paste(samplei, "_seed", seed, "k5", sep="")
        fileName = paste("~/dEploidPaper/validation/errorCount/", prefix, ".withDrop.errorCount", sep ="")
        if ( file.exists(fileName) ){
            tmpInfo = read.table(fileName, header=F)
            for ( i in 1:dim(tmpInfo)[1] ){
                mytable[mytablei,]$prop = tmpInfo$V1[i]
                mytable[mytablei,]$switchError = tmpInfo$V2[i]
                mytable[mytablei,]$missCopy = tmpInfo$V3[i]
                mytable[mytablei,]$dropError = tmpInfo$V4[i]
                mytable[mytablei,]$sampleName = samplei
                mytable[mytablei,]$strainName = paste(samplei,".", i, sep="")
                labStrain = ""
                if ( i == 1 ){
                    labStrain = "hb3"
                    if ( samplei %in% simpleMix) { labStrain = "3d7" }
                } else {
                    labStrain = "7g8"
                    if ( samplei %in% simpleMix) { labStrain = "dd2" }
                }
                mytable[mytablei,]$labStrain = labStrain
                mytablei = mytablei+1
            }
        }
    }
}


deploid2Color= adjustcolor( "red", alpha.f = 0.5)
#deploid3Color = 6
shapeitColor = adjustcolor( "green", alpha.f = 0.5)
beagleColor = adjustcolor( "blue", alpha.f = 0.5)
lassoColor = "purple"
mytable[["color"]] = deploid2Color
#mytable$color[mytable$sampleName %in% c("PG0395-C","PG0396-C","PG0397-C")] = deploid3Color
mytable[["marker"]] = factor(mytable$labStrain)


mytable2 = data.frame ( prop = numeric(0), switchError = numeric(0), missCopy = numeric(0), dropError = numeric(0), sampleName = character(), strainName = character(), labStrain = character(), marker = character(), color = character(), stringsAsFactors = F)
mytablei = 1
print(simpleMix)
for (samplei in sample2){
#    fileName = paste("../../benchMark/beagle/", samplei, ".errorCount", sep ="")
fileName = paste("~/DEploid-Supplementary-Materials/benchMark/beagle/", samplei, ".withDrop.errorCount", sep ="")

    if ( file.exists(fileName) ){
        tmpInfo = read.table(fileName, header=F)
        for ( i in 1:dim(tmpInfo)[1] ){
            mytable2[mytablei,]$prop = tmpInfo$V1[i]
            mytable2[mytablei,]$switchError = tmpInfo$V2[i]
            mytable2[mytablei,]$missCopy = tmpInfo$V3[i]
            mytable2[mytablei,]$dropError = tmpInfo$V4[i]
            mytable2[mytablei,]$sampleName = samplei
            mytable2[mytablei,]$strainName = paste(samplei,".", i, sep="")
            labStrain = ""
            if ( i == 1 ){
                labStrain = "hb3"
                if ( samplei %in% simpleMix) { labStrain = "3d7" }
            } else {
                labStrain = "7g8"
                if ( samplei %in% simpleMix) { labStrain = "dd2" }
            }
            mytable2[mytablei,]$labStrain = labStrain
            mytable2[mytablei,]$color = 5
            mytablei = mytablei+1
        }
    }
}

mytable2[["marker"]] = factor(mytable2$labStrain)


mytable3 = data.frame ( prop = numeric(0), switchError = numeric(0), missCopy = numeric(0), dropError = numeric(0), sampleName = character(), strainName = character(), labStrain = character(), marker = character(), color = character(), stringsAsFactors = F)
mytablei = 1
print(simpleMix)
for (samplei in sample2){
#    fileName = paste("../../benchMark/shapeit/", samplei, ".errorCount", sep ="")
fileName = paste("~/DEploid-Supplementary-Materials/benchMark/shapeit/", samplei, ".withDrop.errorCount", sep ="")
    if ( file.exists(fileName) ){
        tmpInfo = read.table(fileName, header=F)
        for ( i in 1:dim(tmpInfo)[1] ){
            mytable3[mytablei,]$prop = tmpInfo$V1[i]
            mytable3[mytablei,]$switchError = tmpInfo$V2[i]
            mytable3[mytablei,]$missCopy = tmpInfo$V3[i]
            mytable3[mytablei,]$dropError = tmpInfo$V4[i]
            mytable3[mytablei,]$sampleName = samplei
            mytable3[mytablei,]$strainName = paste(samplei,".", i, sep="")
            labStrain = ""
            if ( i == 1 ){
                labStrain = "hb3"
                if ( samplei %in% simpleMix) { labStrain = "3d7" }
            } else {
                labStrain = "7g8"
                if ( samplei %in% simpleMix) { labStrain = "dd2" }
            }
            mytable3[mytablei,]$labStrain = labStrain
            mytable3[mytablei,]$color = 5
            mytablei = mytablei+1
        }
    }
}

mytable3[["marker"]] = factor(mytable3$labStrain)



mytable4 = data.frame ( prop = numeric(0), switchError = numeric(0), missCopy = numeric(0), dropError = numeric(0), sampleName = character(), strainName = character(), labStrain = character(), marker = character(), color = character(), stringsAsFactors = F)
mytablei = 1
print(simpleMix)
for (samplei in sample2){
#    fileName = paste("../../benchMark/shapeit/", samplei, ".errorCount", sep ="")
fileName = paste("~/DEploid-Lasso-prep/testDEploid-Lasso/", samplei, ".withDrop.errorCount", sep ="")
    if ( file.exists(fileName) ){
        tmpInfo = read.table(fileName, header=F)
        for ( i in 1:dim(tmpInfo)[1] ){
            mytable4[mytablei,]$prop = tmpInfo$V1[i]
            mytable4[mytablei,]$switchError = tmpInfo$V2[i]
            mytable4[mytablei,]$missCopy = tmpInfo$V3[i]
            mytable4[mytablei,]$dropError = tmpInfo$V4[i]
            mytable4[mytablei,]$sampleName = samplei
            mytable4[mytablei,]$strainName = paste(samplei,".", i, sep="")
            labStrain = ""
            if ( i == 1 ){
                labStrain = "hb3"
                if ( samplei %in% simpleMix) { labStrain = "3d7" }
            } else {
                labStrain = "7g8"
                if ( samplei %in% simpleMix) { labStrain = "dd2" }
            }
            mytable4[mytablei,]$labStrain = labStrain
            mytable4[mytablei,]$color = 5
            mytablei = mytablei+1
        }
    }
}

mytable4[["marker"]] = factor(mytable4$labStrain)



layout(rbind(1,2,3,4), heights=c(10,10,10,1))  # put legend on bottom 1/8th of the chart
#############################################################################################################################################################################

#par(mfrow = c(2,1))
mytable$switchError[mytable$switchError==0] = 0.05
mytable2$switchError[mytable2$switchError==0] = 0.05
mytable3$switchError[mytable3$switchError==0] = 0.05
mytable4$switchError[mytable4$switchError==0] = 0.05

plot(mytable$prop, mytable$switchError, ylab= "Number of switch errors", xlab="Strain proportions", col=mytable$color, log = "y", ylim= c(0.05, 900), xlim = c(0,1.05), type="n", cex.lab=1.6, cex.axis = 1.4,yaxt="n")
at = c(0.05, 0.09, 0.16, 0.25, 1, 10, 100, 500)
axis(2, at=at, labels=c("0", "", "...", "", "1", "10", "100", "500"), las=0, lwd = 1, cex=2, cex.axis=1.4)
strains = unique(mytable2$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable2$strainName == strain)

    tmpMarker = as.numeric(unique(mytable2$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable2$color[tmpIndex])
    if (mytable2$prop[tmpIndex] > 0 & mytable2$prop[tmpIndex] < 1){
        x = c(x, mean(mytable2$prop[tmpIndex]))
        y = c(y, mean(mytable2$switchError[tmpIndex]))
        points( mean(mytable2$prop[tmpIndex]), mean(mytable2$switchError[tmpIndex]), pch = tmpMarker, col=adjustcolor( beagleColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable2[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.25, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=beagleColor, cex=1.1, lwd=2, lty=1)

strains = unique(mytable3$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable3$strainName == strain)
    tmpMarker = as.numeric(unique(mytable3$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable3$color[tmpIndex])
    if (mytable3$prop[tmpIndex] > 0 & mytable3$prop[tmpIndex] < 1){
        x = c(x, mean(mytable3$prop[tmpIndex]))
        y = c(y, mean(mytable3$switchError[tmpIndex]))
        points( mean(mytable3$prop[tmpIndex]), mean(mytable3$switchError[tmpIndex]), pch = tmpMarker, col=adjustcolor(shapeitColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable3[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.25, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=shapeitColor, cex=1.1, lwd=2, lty=1)

strains = unique(mytable4$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable4$strainName == strain)
    tmpMarker = as.numeric(unique(mytable4$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable4$color[tmpIndex])
    if (mytable4$prop[tmpIndex] > 0 & mytable4$prop[tmpIndex] < 1){
        x = c(x, mean(mytable4$prop[tmpIndex]))
        y = c(y, mean(mytable4$switchError[tmpIndex]))
        points( mean(mytable4$prop[tmpIndex]), mean(mytable4$switchError[tmpIndex]), pch = tmpMarker, col=adjustcolor(lassoColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable4[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.25, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=lassoColor, cex=1.1, lwd=2, lty=1)


strains = unique(mytable$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable$strainName == strain)
    tmpMarker = as.numeric(unique(mytable$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable$color[tmpIndex])
    x = c(x, mean(mytable$prop[tmpIndex]))
    y = c(y, mean(mytable$switchError[tmpIndex]))
    points( mean(mytable$prop[tmpIndex]), mean(mytable$switchError[tmpIndex]), pch = tmpMarker, col=adjustcolor(tmpColor, alpha.f = 0.3), cex=1.1, lwd=2)
}
smoothed = loess.smooth(x,y, span=0.4, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=tmpColor, cex=1.1, lwd=2, lty=1)


legend("topright", legend = c("3D7", "7G8", "Dd2", "HB3"), pch = c(1,0,3,4), cex=1.6, pt.lwd=2)
#legend("top", legend = c("Beagle", "Shapeit", "DEploid (Mixture of 2)", "DEploid (Mixture of 3)"), text.col = c(3,8,2,4), cex=1.4, ncol=2)
legend("top", legend = c("DEploid (Mixture of 2)", "BEAGLE (Mixture of 2)", "SHAPEIT (Mixture of 2)", "DEploid-Lasso"), text.col = c(deploid2Color,beagleColor,shapeitColor, lassoColor), cex=1.4, ncol=2)

#############################################################################################################################################################################


#mytable$missCopy = mytable$missCopy + 100
plot(mytable$prop, mytable$missCopy/totalSites, ylab= "Genotype error rate", xlab="Strain proportions", col=mytable$color, log="y", ylim= c(0.0005, .2), xlim = c(0,1.05), type="n", cex.lab=1.6, cex.axis = 1.4, yaxt = "n")
at = c(0.0005, 5e-3, 5e-2, 0.2, 1)
axis(2, at=at, labels=c("0.0005", "0.005", "0.05", "0.2", "1"), las=0, lwd = 1, cex=2, cex.axis=1.4)
strains = unique(mytable3$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable3$strainName == strain)

    tmpMarker = as.numeric(unique(mytable3$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable3$color[tmpIndex])
    if (mytable3$prop[tmpIndex] > 0 & mytable3$prop[tmpIndex] < 1){
        x = c(x, mean(mytable3$prop[tmpIndex]))
        y = c(y, mean(mytable3$missCopy[tmpIndex])/totalSites)
        points( mean(mytable3$prop[tmpIndex]), mean(mytable3$missCopy[tmpIndex])/totalSites, pch = tmpMarker, col=adjustcolor(shapeitColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable3[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.25, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=shapeitColor, cex=1.1, lwd=2, lty=1)

strains = unique(mytable4$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable4$strainName == strain)

    tmpMarker = as.numeric(unique(mytable4$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable4$color[tmpIndex])
    if (mytable4$prop[tmpIndex] > 0 & mytable4$prop[tmpIndex] < 1){
        x = c(x, mean(mytable4$prop[tmpIndex]))
        y = c(y, mean(mytable4$missCopy[tmpIndex])/18570)
        points( mean(mytable4$prop[tmpIndex]), mean(mytable4$missCopy[tmpIndex])/18570, pch = tmpMarker, col=adjustcolor(lassoColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable4[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.25, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=lassoColor, cex=1.1, lwd=2, lty=1)


strains = unique(mytable2$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable2$strainName == strain)

    tmpMarker = as.numeric(unique(mytable2$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable2$color[tmpIndex])
    if (mytable2$prop[tmpIndex] > 0 & mytable2$prop[tmpIndex] < 1){
        x = c(x, mean(mytable2$prop[tmpIndex]))
        y = c(y, mean(mytable2$missCopy[tmpIndex])/totalSites)
        points( mean(mytable2$prop[tmpIndex]), mean(mytable2$missCopy[tmpIndex])/totalSites, pch = tmpMarker, col=adjustcolor( beagleColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable2[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.25, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=beagleColor, cex=1.1, lwd=2, lty=1)

strains = unique(mytable$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable$strainName == strain)
    tmpMarker = as.numeric(unique(mytable$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable$color[tmpIndex])
#    print(strain)
#    if (as.character(strain) == "PG0402-C.2"){
#        text( mean(mytable$prop[tmpIndex]), mean(mytable$missCopy[tmpIndex])/totalSites, label = "PG0402-C.7G8", col=tmpColor, cex=1.5)
#    }else{
    x = c(x, mean(mytable$prop[tmpIndex]))
    y = c(y, mean(mytable$missCopy[tmpIndex])/totalSites)
        points( mean(mytable$prop[tmpIndex]), mean(mytable$missCopy[tmpIndex])/totalSites, pch = tmpMarker, col=adjustcolor(tmpColor, alpha.f = 0.3), cex=1.1, lwd=2)
#    }
}
smoothed = loess.smooth(x,y, span=0.21, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=tmpColor, cex=1.1, lwd=2, lty=1)
legend("topright", legend = c("3D7", "7G8", "Dd2", "HB3"), pch = c(1,0,3,4), cex=1.6, pt.lwd=2)
#legend("top", legend = c("Mixture of 2", "Mixture of 3", "Beagle", "Shapeit"), text.col = c(2,4,3,8), cex=1.4, ncol=4)
legend("top", legend = c("DEploid (Mixture of 2)", "BEAGLE (Mixture of 2)", "SHAPEIT (Mixture of 2)", "DEploid-Lasso"), text.col = c(deploid2Color,beagleColor,shapeitColor, lassoColor), cex=1.4, ncol=2)

#############################################################################################################################################################################

mytable$dropError[mytable$dropError==0] = 1
mytable2$dropError[mytable2$dropError==0] = 1
mytable3$dropError[mytable3$dropError==0] = 1
mytable4$dropError[mytable4$dropError==0] = 1

plot(mytable$prop, mytable$dropError, ylab= "Drop-out error fraction", xlab="Strain proportions", col=mytable$color, log = "y", ylim= c(5.704507e-05, 2), xlim = c(0,1.05), type="n", cex.lab=1.6, cex.axis = 1.4, yaxt = "n")
at = c(5.704507e-05, 0.000100105, 0.000150105, 0.000200105, 1e-3, 1e-2, 1e-1, 1)
axis(2, at=at, labels=c("0", "", "...", "", "0.001", "0.01", "0.1", "1"), las=0, lwd = 1, cex=2, cex.axis=1.4)


strains = unique(mytable2$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable2$strainName == strain)

    tmpMarker = as.numeric(unique(mytable2$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable2$color[tmpIndex])
    if (mytable2$prop[tmpIndex] > 0 & mytable2$prop[tmpIndex] < 1){
        x = c(x, mean(mytable2$prop[tmpIndex]))
        y = c(y, mean(mytable2$dropError[tmpIndex])/totalSites)
        points( mean(mytable2$prop[tmpIndex]), mean(mytable2$dropError[tmpIndex])/totalSites, pch = tmpMarker, col=adjustcolor( beagleColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable2[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.25, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=beagleColor, cex=1.1, lwd=2, lty=1)



strains = unique(mytable3$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable3$strainName == strain)

    tmpMarker = as.numeric(unique(mytable3$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable3$color[tmpIndex])
    if (mytable3$prop[tmpIndex] > 0 & mytable3$prop[tmpIndex] < 1 ){
        x = c(x, mean(mytable3$prop[tmpIndex]))
        y = c(y, mean(mytable3$dropError[tmpIndex])/totalSites)
        points( mean(mytable3$prop[tmpIndex]), mean(mytable3$dropError[tmpIndex])/totalSites, pch = tmpMarker, col=adjustcolor(shapeitColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable3[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.21, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=shapeitColor, cex=1.1, lwd=2, lty=1)


strains = unique(mytable4$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable4$strainName == strain)

    tmpMarker = as.numeric(unique(mytable4$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable4$color[tmpIndex])
    if (mytable4$prop[tmpIndex] > 0 & mytable4$prop[tmpIndex] < 1 ){
        x = c(x, mean(mytable4$prop[tmpIndex]))
        y = c(y, mean(mytable4$dropError[tmpIndex])/totalSites)
        points( mean(mytable4$prop[tmpIndex]), mean(mytable4$dropError[tmpIndex])/totalSites, pch = tmpMarker, col=adjustcolor(lassoColor, alpha.f = 0.3), cex=1.1, lwd=2)
    } else {
        print(mytable4[tmpIndex,])
    }
}
smoothed = loess.smooth(x,y, span=0.21, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=lassoColor, cex=1.1, lwd=2, lty=1)


strains = unique(mytable$strainName)
x = c()
y = c()
for ( strain in strains ){
    tmpIndex = which(mytable$strainName == strain)
    tmpMarker = as.numeric(unique(mytable$marker[tmpIndex]))
    if ( tmpMarker == 2) {tmpMarker = 0}
    tmpColor = unique(mytable$color[tmpIndex])
    x = c(x, mean(mytable$prop[tmpIndex]))
    y = c(y, mean(mytable$dropError[tmpIndex])/totalSites)
    points( mean(mytable$prop[tmpIndex]), mean(mytable$dropError[tmpIndex])/totalSites, pch = tmpMarker, col=adjustcolor(tmpColor, alpha.f = 0.3), cex=1.1, lwd=2)
}
smoothed = loess.smooth(x,y, span=0.21, family = c("gaussian"))
lines(smoothed$x, smoothed$y, col=tmpColor, cex=1.1, lwd=2, lty=1)


legend("topright", legend = c("3D7", "7G8", "Dd2", "HB3"), pch = c(1,0,3,4), cex=1.6, pt.lwd=2)
#legend("top", legend = c("Beagle", "Shapeit", "DEploid (Mixture of 2)", "DEploid (Mixture of 3)"), text.col = c(3,8,2,4), cex=1.4, ncol=2)
legend("top", legend = c("DEploid (Mixture of 2)", "BEAGLE (Mixture of 2)", "SHAPEIT (Mixture of 2)", "DEploid-Lasso"), text.col = c(deploid2Color,beagleColor,shapeitColor, lassoColor), cex=1.4, ncol=2)

dev.off()


library(dplyr)

mytable_grouped = group_by(mytable, strainName) %>% summarise(., perfect.panel.switchError.m = mean(switchError),
                                               perfect.panel.missCopy.m = mean(missCopy),
                                               perfect.panel.dropError.m = mean(dropError) )

merged.tab = merge(mytable_grouped, mytable4)

for (i in c("switchError", "missCopy", "dropError")){
    x = merged.tab[[i]]
    y = merged.tab[[paste("perfect.panel.", i, ".m", sep = "")]]
    png(paste("pefectPanelVsDEploidLasso.", i, ".png", sep = ""))
    par(mar=c(5,5,2,2))
    plot(x, y, log="xy", xlim = range(x,y), ylim=range(x,y), xlab = "Perfect panel", ylab = "DEploid-LASSO", main = paste(i, "corr:", round(cor(x,y), digits = 3)), cex.lab = 1.5, cex.main = 1.5)
    lm1 = lm(y~x)
    newx <- seq(min(x), max(x), length.out=100)
    preds <- predict(lm1, newdata = data.frame(x=newx))
    dev.off()
#    lines(newx, preds)
}

