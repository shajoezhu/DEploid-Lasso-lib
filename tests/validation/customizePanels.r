rm(list=ls())
library(DEploid)
library(dplyr)

common_plot <- function(){
    plot(myfit$df, main = "df", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5, col = points_col)
    abline(v = dev.ratio.rate.cut, col = "blue" )

    plot.glmnet(myfit, xvar="norm", col = useColor, main = "norm", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5)
    plot(myfit$a0, main = "a0", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5, col = points_col)
    abline(v = dev.ratio.rate.cut, col = "blue" )

    plot.glmnet(myfit, xvar="lambda", col = useColor, main = "lambda", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5)
    plot(myfit$lambda, main = "lambda", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5, col = points_col)
    abline(v = dev.ratio.rate.cut, col = "blue" )

    plot.glmnet(myfit, xvar="dev", col = useColor, main = "dev.ratio", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5, xlim = c(0,1))
    plot(myfit$dev.ratio, main = "dev.ratio", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5, col = points_col, ylim=c(0,1))
    abline(v = dev.ratio.rate.cut, col = "blue" )

    plot(c(0,0), c(1,1), type = "n", xlab = "", ylab = "", axes=F, main = "", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5)
    plot(c(1:n_lambda), c(dev.ratio.rate,0), main = "dev.ratio change", cex.main = 3, cex.lab = 1.5, cex.axis = 1.5, ylim=c(0,0.02), col = points_col)
    abline(v = dev.ratio.rate.cut, col = "blue" )
}

panel = read.table("panel_from_crosses.txt", header=T, comment.char="", check.names=F)
panelCol = names(panel)
panel = panel[, -which(panelCol %in% c("PG0051-C", "PG0052-C", "7G8", "GB4", "PG0004-CW", "PG0008-CW"))]

strainColors = rep("grey", length(names(panel)))
strainColors[names(panel) %in% read.table("cross_3D7xHB3", header=F)$V1] = "red"
strainColors[names(panel) %in% read.table("cross_7G8xGB4", header=F)$V1] = "blue"
strainColors[names(panel) %in% read.table("cross_HB3xDD2", header=F)$V1] = "green"

i_arr = c(389:415)
correct_df = c(rep(1,6),
               rep(3,3),
               1,
               rep(2,16),
               1)

#n_lambda = 100
dev.explained = c()
for ( idx_i in 1:27 ){
    i = i_arr[idx_i]
    sample_i = paste("PG0", i, "-C", sep="")
#    sample_i = "PG0396-C"
    coverage = extractCoverageFromVcf(paste("jasonMixed/", sample_i, ".vcf.gz", sep = ""))
dev.explained.sample = c()
    for ( chrom in unique(coverage$CHROM) ){
        chrom_idx = which(coverage$CHROM == chrom)
        panelName = paste("crossProgenyOnlyPanel_", chrom, ".txt", sep="")
        if (idx_i == 1){
            write.table(panel[chrom_idx,], file = panelName, sep = "\t", row.names=F, col.names=F)
        }

        wsaf = computeObsWSAF(alt = coverage$altCount[chrom_idx], ref = coverage$refCount[chrom_idx])
        wsafName = paste("tmp/", sample_i, "_", chrom, ".wsaf", sep="")
        write.table(wsaf, file = wsafName, row.names=F, col.names=F)

        outFile = paste("tmp/", sample_i, "_", chrom, ".lasso", sep="")
        system (paste("lasso", panelName, wsafName, ">", outFile))

#        inFileName = args[1]

        system(paste("grep TABLE ", outFile, " > ", outFile, ".TABLE", sep=""))
        system(paste("grep BETA ", outFile, " > ", outFile, ".BETA", sep=""))

        mytable = read.table(paste(outFile, ".TABLE", sep=""), header=T, stringsAsFactors = F)[,-1]
        mybeta = read.table(paste(outFile, ".BETA", sep=""), header=T, stringsAsFactors = F)[,-1]

        dev.ratio.rate = diff(mytable$rsq) / mytable$rsq[-1]
        dev.ratio.rate.cut = min(which(dev.ratio.rate<0.001)[1], length(dev.ratio.rate), na.rm = T)

        png( paste(sample_i, "_", chrom, ".png", sep=""), width = 1000, height = 1000)
        par(mar=c(5,5,2,2))
        par(mfrow = c(4,1))
        plot(mytable$df, ylab = "# of predictor variable", cex.lab = 2, cex.axis = 2)
        abline(v = dev.ratio.rate.cut, col = "blue" )

        plot(mytable$rsq, ylab = "Deviance ratio", cex.lab = 2, cex.axis = 2)
        abline(v = dev.ratio.rate.cut, col = "blue" )

        plot(diff(mytable$rsq), ylab = "Deviance ratio differences", cex.lab = 2, cex.axis = 2 )
        abline(v = dev.ratio.rate.cut, col = "blue" )

        v_idx = c()
        plot(c(0, max(mytable$L1norm)),
             c(min(mybeta[,-1]), max(mybeta[,-1])), type="n", xlab = "L1 Norm", ylab = "Coefficient", cex.lab = 2, cex.axis = 2)
        for (i in 1:dim(mybeta)[1]){
             lines(c(0, mytable$L1norm), c(0,mybeta[i,-1]), col = strainColors[mybeta[i,1]])
        }
        abline(v = mytable$L1norm[dev.ratio.rate.cut], col = "blue" )

        legend("left", legend=c("3D7xHB3", "7G8xGB4", "HB3xDD2"), col = c("red", "blue", "green"), lty=1, bty="n", cex = 2)
        dev.off()

        CHROM = coverage$CHROM[chrom_idx]
        POS = coverage$POS[chrom_idx]
        new_tmp_panel = cbind(CHROM, POS, panel[chrom_idx,mybeta[which(mybeta[,1+dev.ratio.rate.cut] > 0),1]] )
        write.table(new_tmp_panel, file = paste("panels/", sample_i, "_", chrom, ".panel.txt", sep =""), quote = F, row.names=F, sep = "\t")

    }
}


