####
## Script to preform permutation test on the specified matrices of gene CN counts
####

setwd("/data/proj/chilense/30_genomes_outputs/kai_wei/variant_call/CNVs/metaSV/control-freec/CN/")

#normalization of read depth by median value

data <- read.table("gene_read_depth.txt", header = T, row.names = 1)

x=data.frame(data)

nor <- apply(x,2,function(x){(x/median(x))*2})


write.table(nor, file = "gene_read_depth_normolized.txt")


library(devtools)
library(Hmisc)


cnvCN = read.table(file="gene_CN.txt", row.names = 1, header = T)
cnvRD = read.table(file="gene_read_depth_normolized.txt", row.names = 1, header = T)


#define sample lists

samp=c("LA1963_t1","LA1963_t3","LA1963_t5","LA1963_t7","LA1963_t9","LA2931_t2","LA2931_t3","LA2931_t4","LA2931_t5","LA2931_t6",
       "LA2932_1","LA2932_8","LA2932_12","LA2932_20","LA2932_22","LA3111_t3","LA3111_t5","LA3111_t9","LA3111_t10","LA3111_t15",
       "LA4107_3","LA4107_6","LA4107_9","LA4107_t5","LA4107_t11","LA4117_1","LA4117_4","LA4117_5","LA4117_10","LA4117_t15",
       "LA4330_t1","LA4330_t4","LA4330_t6","LA4330_t9","LA4330_t12")

spec=c("LA1963","LA1963","LA1963","LA1963","LA1963","LA2931","LA2931","LA2931","LA2931","LA2931","LA2932","LA2932","LA2932","LA2932",
       "LA2932","LA3111","LA3111","LA3111","LA3111","LA3111","LA4107","LA4107","LA4107","LA4107","LA4107","LA4117","LA4117","LA4117",
       "LA4117","LA4117","LA4330","LA4330","LA4330","LA4330","LA4330")


#calculate per-gene CN variance in BB, PB and BBPB
VCNall = t(apply(cnvCN, 1, function(y) c(var(y[1:5]), var(y[6:10]), var(y[11:15]), var(y[16:20]), var(y[21:25]), var(y[26:30]),
                                         var(y[31:35]), var(y[1:35])))) ##equal var.p in Excel

VRDall = t(apply(cnvRD, 1, function(y) c(var(y[1:5]), var(y[6:10]), var(y[11:15]), var(y[16:20]), var(y[21:25]), var(y[26:30]),
                                         var(y[31:35]), var(y[1:35])))) ##equal var.p in Excel

#calculate per-gene Vst
VstCNall = (VCNall[,8] - ((VCNall[,1]*5 + VCNall[,2]*5 + VCNall[,3]*5 + VCNall[,4]*5 + VCNall[,5]*5 + VCNall[,6]*5 + VCNall[,7]*5)/35))/VCNall[,8]

VstRDall = (VRDall[,8] - ((VRDall[,1]*5 + VRDall[,2]*5 + VRDall[,3]*5 + VRDall[,4]*5 + VRDall[,5]*5 + VRDall[,6]*5 + VRDall[,7]*5)/35))/VRDall[,8]

###
#Remove all genes with Vst <= 0 to create a reduced matric containing on CN variable genes
###

NonzeroVstCN = names(which(VstCNall>0))
NonzeroVstRD = names(which(VstRDall>0))

cnvCNreduced = cnvCN[(rownames(cnvCN) %in% NonzeroVarCN),]
cnvRDreduced = cnvRD[(rownames(cnvRD) %in% NonzeroVarRD),]

####
#calculate Vst (this is redundant from above but it only takes a few seconds)

V1 = t(apply(cnvCNreduced, 1, function(y) c(var(y[1:5]), var(y[6:10]), var(y[11:15]), var(y[16:20]), var(y[21:25]), var(y[26:30]),
                                     var(y[31:35]), var(y[1:35]))))

V2 = t(apply(cnvRDreduced, 1, function(y) c(var(y[1:5]), var(y[6:10]), var(y[11:15]), var(y[16:20]), var(y[21:25]), var(y[26:30]),
                                            var(y[31:35]), var(y[1:35]))))

VstCNreduced = (V1[,8] - ((V1[,1]*5 + V1[,2]*5 + V1[,3]*5 + V1[,4]*5 + V1[,5]*5 + V1[,6]*5 + V1[,7]*5)/35))/V1[,8]

VstRDreduced = (V2[,8] - ((V2[,1]*5 + V2[,2]*5 + V2[,3]*5 + V2[,4]*5 + V2[,5]*5 + V2[,6]*5 + V2[,7]*5)/35))/V2[,8]

#shuffle reduced CN matrix
cnvCNrdm = cnvCNreduced[, c(1:2, sample(3:ncol(cnvCNreduced)))]
cnvRDrdm = cnvRDreduced[, c(1:2, sample(3:ncol(cnvRDreduced)))]

#calculate Vst on shuffled matrix

V1rdm = t(apply(cnvCNrdm, 1, function(y) c(var(y[1:5]), var(y[6:10]), var(y[11:15]), var(y[16:20]), var(y[21:25]), var(y[26:30]),
                                        var(y[31:35]), var(y[1:35]))))

V2rdm = t(apply(cnvRDrdm, 1, function(y) c(var(y[1:5]), var(y[6:10]), var(y[11:15]), var(y[16:20]), var(y[21:25]), var(y[26:30]),
                                           var(y[31:35]), var(y[1:35]))))

VstCNrdm = (V1rdm[,8] - ((V1rdm[,1]*5 + V1rdm[,2]*5 + V1rdm[,3]*5 + V1rdm[,4]*5 + V1rdm[,5]*5 + V1rdm[,6]*5 + V1rdm[,7]*5)/35))/V1rdm[,8]
VstRDrdm = (V2rdm[,8] - ((V2rdm[,1]*5 + V2rdm[,2]*5 + V2rdm[,3]*5 + V2rdm[,4]*5 + V2rdm[,5]*5 + V2rdm[,6]*5 + V2rdm[,7]*5)/35))/V2rdm[,8]

VstPerCN = VstCNrdm
VstPerRD = VstRDrdm


#repeat x1000
for (i in 1:999){
  
  cnvCNrdm = cnvCNreduced[, c(1:2, sample(3:ncol(cnvCNreduced)))]
  cnvRDrdm = cnvRDreduced[, c(1:2, sample(3:ncol(cnvRDreduced)))]
  
  V1rdm = t(apply(cnvCNrdm, 1, function(y) c(var(y[1:5]), var(y[6:10]), var(y[11:15]), var(y[16:20]), var(y[21:25]), var(y[26:30]),
                                             var(y[31:35]), var(y[1:35]))))
  
  V2rdm = t(apply(cnvRDrdm, 1, function(y) c(var(y[1:5]), var(y[6:10]), var(y[11:15]), var(y[16:20]), var(y[21:25]), var(y[26:30]),
                                             var(y[31:35]), var(y[1:35]))))
  
  VstCNrdm = (V1rdm[,8] - ((V1rdm[,1]*5 + V1rdm[,2]*5 + V1rdm[,3]*5 + V1rdm[,4]*5 + V1rdm[,5]*5 + V1rdm[,6]*5 + V1rdm[,7]*5)/35))/V1rdm[,8]
  VstRDrdm = (V2rdm[,8] - ((V2rdm[,1]*5 + V2rdm[,2]*5 + V2rdm[,3]*5 + V2rdm[,4]*5 + V2rdm[,5]*5 + V2rdm[,6]*5 + V2rdm[,7]*5)/35))/V2rdm[,8]
  
  #add vector to permutation test matrix
  VstPerCN = cbind(VstCNrdm, VstPerCN)
  VstPerRD = cbind(VstRDrdm, VstPerRD)
}

############
##Extract cutoffs for each gene from permutations
############
VstCNcutoff.95 = apply(VstPerCN, 1, function(x) quantile(x, probs=.95))
VstCNcutoff.99 = apply(VstPerCN, 1, function(x) quantile(x, probs=.99))

VstRDcutoff.95 = apply(VstPerRD, 1, function(x) quantile(x, probs=.95))
VstRDcutoff.99 = apply(VstPerRD, 1, function(x) quantile(x, probs=.99))

VstCutoffsCN = cbind( VstCNcutoff.95, VstCNcutoff.99, VstCNreduced)
VstCutoffsRD = cbind( VstRDcutoff.95, VstRDcutoff.99, VstRDreduced)







