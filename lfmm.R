#devtools::install_github("bcm-uga/lfmm")
#devtools::install_github("bcm-uga/LEA")
#devtools::install_github("bcm-uga/TESS3_encho_sen")
setwd("/data/proj/chilense/30_genomes_outputs/kai_wei/variant_call/CNVs/metaSV/VST/LFMM")
library(lfmm)
library(LEA)
library(qvalue)
library(tess3r)
library(RSpectra)

cn <- read.table('CN_Cutoff95.txt', header = T, row.names = 1)

# a principal component analysis (PCA) can reveal some ‘structure’ in the copy number data
pc <- prcomp(cn)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(4,pc$sdev[4]^2, type = "h", lwd = 3, col = "blue") ## there are around K=4 main components in the data, same with admixture

clim <- read.table("climate.txt", header = T)
Y <- scale(cn)
X <- scale(as.numeric(clim$PETColdestQuarter))


axex.lfmm <- lfmm_ridge(Y = Y,
                        X = X, 
                        K = 4)


PETColdestQuarter <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = axex.lfmm, 
                calibrate = "gif")

pvalues <- pv$calibrated.pvalue

par(mfrow=c(5,1),mar=c(4,4,1,1))
hist(PETDriestQuarter$pvalue, breaks = 20, main = "PETDriestQuarter", xlab = "p-values")
qqplot(rexp(length(PETDriestQuarter$calibrated.pvalue), rate = log(10)),
       -log10(PETDriestQuarter$calibrated.pvalue), ylab = "Observed -log10 (P)", xlab="Expected -log10 (P)",
       pch = 19, cex = .4, col="red", main="PETDriestQuarter")

abline(0,1)


L = ncol(cn)

plot(-log10(PETDriestQuarter$calibrated.pvalue), 
     pch = 19, 
     cex = .3,
     xlab = "Genes", ylab = "-log P",
     col = "grey")


write.table(fdr, file="test.txt")

fdr <- p.adjust(bio7$calibrated.pvalue, method = "BH", n = length(bio7$calibrated.pvalue))

fdr <- round(p.adjust(bio7$calibrated.pvalue, "BH"), 5)
