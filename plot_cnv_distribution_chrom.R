setwd("/data/proj/chilense/30_genomes_outputs/kai_wei/variant_call/CNVs/merge/merged_samples/")

library(karyoploteR)
del <- read.table("merged_DEL.txt", header = T)
dup <- read.table("merged_DUP.txt", header = T)

par(omi=c(0,0.5,0,0))

custom.genome <- toGRanges("chrom.txt")

##plot DUP
kp <- plotKaryotype(genome = custom.genome, plot.type=4)
grgenerange = with(dup, GRanges(Chrom, IRanges(start=Start, end=End)))
kpPlotDensity(kp, grgenerange, window.size = 2e6, col="#00BFC4")
kpAxis(kp, ymin = 0, ymax=330)


##plot DEL
kp <- plotKaryotype(genome = custom.genome, plot.type=4)
grgenerange = with(del, GRanges(Chrom, IRanges(start=Start, end=End)))
kpPlotDensity(kp, grgenerange, window.size = 2e6, col="#F8766D")
kpAxis(kp, ymin = 0, ymax=660)



