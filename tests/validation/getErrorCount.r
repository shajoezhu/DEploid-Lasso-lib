#!/usr/bin/env Rscript
rm(list=ls())

args = (commandArgs(TRUE))

prefix = args[1]

source("common.r")

panel = read.table("labStrains.eg.panel.txt", header=T)
#keepingIdx = read.table("keeping.idx.txt", header= T)

#panel = panel[keepingIdx$IDX,]
#Need to read exclude sites ... first trim the panel!

#panel = read.table("labStrains.eg.panel.txt", header=T)

endAt = cumsum(table(panel[,1]))
beginAt = c(1, 1+endAt[-length(endAt)])

chromLength = (endAt - beginAt+1)

Ref1 = panel[,5] # HB3
Ref1Name = "HB3"
Ref2 = panel[,6] # 7G8
Ref2Name = "7G8"

panelID = ""
#panelID = ".asiaAfrica"

true.prop = list()
# 3d7, dd2, 7g8, hb3
true.prop[["PG0389-C"]] = c(0.9, 0.1, 0, 0)
true.prop[["PG0390-C"]] = c(0.8, 0.2, 0, 0)
true.prop[["PG0391-C"]] = c(0.67, 0.33, 0, 0)
true.prop[["PG0392-C"]] = c(0.33, 0.67, 0, 0)
true.prop[["PG0393-C"]] = c(0.2, 0.8,0,0)
true.prop[["PG0394-C"]] = c(0.1, 0.9,0,0)
true.prop[["PG0395-C"]] = c(0, 0.33, 0.34, 0.33)
true.prop[["PG0396-C"]] = c(0, 0.25, 0.5, 0.25)
true.prop[["PG0397-C"]] = c(0, 0.143, 0.714, 0.143)
true.prop[["PG0398-C"]] = c(0, 0, 0, 1)
true.prop[["PG0399-C"]] = c(0, 0, 0.01, 0.99)
true.prop[["PG0400-C"]] = c(0, 0, 0.05, 0.95)
true.prop[["PG0401-C"]] = c(0, 0, 0.1, 0.9)
true.prop[["PG0402-C"]] = c(0, 0, 0.15, 0.85)
true.prop[["PG0403-C"]] = c(0, 0, 0.2, 0.8)
true.prop[["PG0404-C"]] = c(0, 0, 0.25, 0.75)
true.prop[["PG0405-C"]] = c(0, 0, 0.3, 0.7)
true.prop[["PG0406-C"]] = c(0, 0, 0.4, 0.6)
true.prop[["PG0407-C"]] = c(0, 0, 0.5, 0.5)
true.prop[["PG0408-C"]] = c(0, 0, 0.6, 0.4)
true.prop[["PG0409-C"]] = c(0, 0, 0.7, 0.3)
true.prop[["PG0410-C"]] = c(0, 0, 0.75, 0.25)
true.prop[["PG0411-C"]] = c(0, 0, 0.80, 0.2)
true.prop[["PG0412-C"]] = c(0, 0, 0.85, 0.15)
true.prop[["PG0413-C"]] = c(0, 0, 0.95, 0.05)
true.prop[["PG0414-C"]] = c(0, 0, 0.99, 0.01)
true.prop[["PG0415-C"]] = c(0, 0, 1,0)

i_arr = c(389:394, 399:414)

for ( idx_i in 1:22 ){
    i = i_arr[idx_i]
    sample_i = paste("PG0", i, "-C", sep="")
#    prop.tmp = paste(true.prop[[sample_i]][true.prop[[sample_i]]>0], collapse = " ")

#sampleName = "PG0393-C"
prefix = sample_i
#for ( seed in 1:15 ){

#    prefix = paste("PG0390-C/PG0390-C_seed", seed, "k2", sep="")
    if ( sample_i %in% c("PG0389-C", "PG0390-C", "PG0391-C", "PG0392-C", "PG0393-C", "PG0394-C") ){
        Ref1 = panel[,3] # 3d7
        Ref1Name = "3d7"
        Ref2 = panel[,4] # Dd2
        Ref2Name = "Dd2"
    } else {
        Ref1 = panel[,5] # HB3
        Ref1Name = "HB3"
        Ref2 = panel[,6] # 7G8
        Ref2Name = "7G8"
    }

    png(paste(prefix, "compareHap.png", sep=""), width = 1920, height = 1080)
    ncol = ceiling(length(endAt)/2)
    par(mfrow = c(ncol,length(endAt)/ncol))
#    tmpProp = read.table(paste(prefix,".prop",sep=""), header=F)
#    prop = as.numeric(tmpProp[dim(tmpProp)[1],])
prop = true.prop[[sample_i]][true.prop[[sample_i]]>0]
    hap = as.matrix(read.table(paste("haps/", prefix, panelID, ".hap",sep=""), header=T)[,c(-1,-2)])

#    colIndex = which(prop>0.01)
    colIndex = c(1,2)

    prop.corrected = prop[colIndex]
    hap.corrected = hap#[,colIndex,drop=FALSE]

    printed.prop = prop.corrected[getIndex2(hap.corrected, Ref1, Ref2)]

    switchError = c(0, 0)
    mutError = c(0, 0)
    dropError = c(0, 0)
    for ( chrom in 1:length(beginAt)){
        tmpHap = hap.corrected[beginAt[chrom]:endAt[chrom],,drop=FALSE]
#        tmpProp = prop.corrected
        tmpRef1 = Ref1[beginAt[chrom]:endAt[chrom]]
        tmpRef2 = Ref2[beginAt[chrom]:endAt[chrom]]

        if ( length(prop.corrected) == 2 ){
            rearranged.Index = getIndex2(tmpHap, tmpRef1, tmpRef2)
            tmpHap = tmpHap[,rearranged.Index,drop=FALSE]
#            tmpProp = prop.corrected[rearranged.Index]
        }

        tmpRef = cbind(tmpRef1, tmpRef2)
        hapAndError = measure.error.joe.with.drop.2strain(t(tmpHap), t(tmpRef),2)

        tmpTitle = paste(rownames(table(panel[,1]))[chrom], sum(hapAndError$switchError), "switch errors", sum(hapAndError$mutError), "miss copy errors")

#        fun.plotHapWithProp (hapAndError$hap, prop,
#             tmpTitle,
#             max(chromLength))
	cat("switchError ", switchError, "\n")
	cat("hapAndError$switchError ", hapAndError$switchError,"\n")
        switchError = switchError + hapAndError$switch.error
	cat("mutError ", mutError,"\n")
	cat("hapAndError$mutError ", hapAndError$mutError,"\n")
        mutError = mutError + hapAndError$genotype.error
        dropError = dropError + hapAndError$drop.error

    }
#    if ( length(prop.corrected) == 2 ){
        printed.prop = sort(printed.prop, decreasing=T)[sort.int(mutError+dropError, decreasing=F, index.return=T)$ix]
        write.table(cbind(printed.prop, switchError, mutError, dropError), file = paste(prefix,panelID,".withDrop.errorCount", sep=""), quote = F, row.names=F, col.names=F)
#    }
    dev.off()

}
warnings()
