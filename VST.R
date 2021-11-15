setwd("/data/proj/chilense/30_genomes_outputs/kai_wei/variant_call/CNVs/metaSV/VST")

data <- read.table("LA2932_LA4107.txt", header=T)


vst <- (apply(data,1,var)- (apply(data[,1:5],1,var)*5 + apply(data[,6:10],1,var)*5)/10)/apply(data,1,var)

vst[vst<0] <- 0
mean(vst,  na.rm = T)
sd(vst,  na.rm = T)

write.table(vst, file = "VST_LA2932_LA4107.txt")







