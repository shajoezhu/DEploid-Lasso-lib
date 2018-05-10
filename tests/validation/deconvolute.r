rm(list=ls())
library(DEploid)

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


PLAF = read.table("labStrains.eg.PLAF.txt", header=T, stringsAsFactors=F)

chroms = unique(PLAF$CHROM)
CHROM = PLAF$CHROM
POS = PLAF$POS

i_arr = c(389:415)
panel = ""
#panel = ".asiaAfrica"
for ( idx_i in 19:19 ){
    i = i_arr[idx_i]
    sample_i = paste("PG0", i, "-C", sep="")
    prop.tmp = paste(true.prop[[sample_i]][true.prop[[sample_i]]>0], collapse = " ")
    haps = c()
    for ( chrom in chroms ){
        cmd = paste("-ref", paste("coverage/", sample_i, "_", chrom, ".ref", sep =""),
                    "-alt", paste("coverage/", sample_i, "_", chrom, ".alt", sep =""),
                    "-plaf", paste("plafs/", chrom, ".plaf", sep = ""),
                    "-panel", paste("panels/", sample_i, "_", chrom, panel, ".panel.txt", sep =""),
                    "-initialP", prop.tmp)
        print(cmd)
        dEploidout = dEploid(cmd)
        haps = rbind(haps, t(dEploidout$Haps))
    }
    colnames(haps) = paste("h", 1:length(true.prop[[sample_i]][true.prop[[sample_i]]>0]), sep="")
    write.table(cbind(CHROM, POS, haps), file = paste("haps/", sample_i, panel, ".hap", sep=""), row.names=F, quote=F, sep="\t")
}
