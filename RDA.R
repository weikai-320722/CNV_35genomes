
# the vegan RDA analysis was performed using gene copy number and climatic data

setwd('//data/proj/chilense/30_genomes_outputs/kai_wei/variant_call/CNVs/metaSV/VST/RDA')
library(vegan)
library(PerformanceAnalytics)

## read copy number and climatic variables
climate <-read.table("climatic_variables.txt",header = T) # read variables including populations and individuals
env1 <- climate[,3:38] # extract variables excluding populations and individuals
gen <- read.table("RDA_input_CN_Cutoff95.txt",row.names = 1, header=T) # read copy number 


## run RDA analysis

# run RDA 
chilense1.rda <- rda(gen ~1, data= env1, scale = T)#one by one analysis between genes copy number and climatic variables
chilense2.rda <- rda(gen ~., data= env1, scale = T)#RDA analysis
chilense5.rda <- ordiR2step (chilense1.rda, scope = formula (chilense2.rda))#Choose a Model by Permutation Tests in Constrained Ordination
vif.cca(chilense5.rda)# analyse linear dependencies among constraints and conditions, remove values over 10 indicating redundant constraints and
summary(eigenvals(chilense5.rda))  # check proportion Explained in every Eigenvalue (RDA)

#  calculate the adjusted R2 to estimate how many variantion can be explained by constrained ordination
RsquareAdj(chilense.rda)

screeplot(chilense5.rda) #plot each constrained axes (RDA) can explained variations

## plot the RDA

levels(climate$Population) <- c("LA1963", "LA2931", "LA2932", "LA3111", "LA4107","LA4117A", "LA4330") # grouping populations
pop <- climate$Population
bg <- c("#F8766D", "#B79F00","#00BA38","#619CFF","#0000FF","#CC0033","#F564E3")# set up color for each population
plot(chilense5.rda, type="n", choices=c(1,2), scaling=1) # plot rda1 and rda2
points(chilense5.rda, display = "species", pch=20, cex=0.7, col="gray", scaling=1) # show genes
points(chilense5.rda, display = "sites", pch=21, cex=1.3, col= NULL,  bg=bg[pop], scaling=1) #show populations
text(chilense5.rda,  display = "bp", col="black", cex=0.8, scaling=1) # show climatic variables
legend("bottomright", legend = levels(pop), bty="n", col="gray", pch=21, cex=1, pt.bg=bg)


## extract outlier genes indicating high correlation between copy number and climatic variables for RDA1 and RDA2

load.rda <- scores(chilense5.rda, choices = c(1:2), display = "species") #extract loading value from RDA! and RDA2
hist(load.rda[,1],xlab = "load.RDA1", main = NULL) # plot loading distribution of constrained axes
hist(load.rda[,2],xlab = "load.RDA2", main = NULL)

# define the function here as outliers, where x is the vector of loadings and z is the number of standard deviations to use
outliers <- function(x,z){
  lims<- mean(x) + c(-1,1) * z* sd(x) # find loadings +/- sd from mean loading
  x[x < lims[1] | x > lims[2]] # locus names in these tail s
}

#  extract significant snps form rda1 and rda2 
cand1 <- outliers(load.rda[,1],2.5) # p-value is about 0.012
cand2 <- outliers(load.rda[,2],2.5)

ncand <- length(cand1) + length(cand2) # summary counts of candidate snps
ncand


#organize our results by making one data frame with the axis, gene name, loading, & correlation with each predictor
cand1 <- cbind.data.frame(rep(1, times=length(cand1)), name=(cand1), unname=(cand1))
cand2 <- cbind.data.frame(rep(2, times=length(cand2)), name=(cand2), unname=(cand2))

colnames(cand1)<-colnames(cand2)<-c("axis", "genes", "loading")
cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)


#Let’s add in the correlations of each candidate genes with the eight environmental predictors:
foo <- matrix(nrow = (ncand), ncol = 5) # ncol is number of variables
colnames(foo) <- c("PETColdestQuarter","Bio8","Bio7","annualPET","PETDriestQuarter")

verna <- subset(env1, select = c("PETColdestQuarter","Bio8","Bio7","annualPET","PETDriestQuarter"))

for (i in 1:length(cand$genes)) {
  nam <- cand[i,1] # nam is snp ID
  snp.gen <- gen[,nam]
  foo[i,] <- apply(verna,2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo) # combine candidates and correlations 

#see which of the variables each candidate SNP is most strongly correlated with
for (i in 1:length(cand$genes)) {
  bar <- cand[i,]
  cand[i,9] <- names(which.max(abs(bar[4:8]))) # gives the variable
  cand[i,10] <- max(abs(bar[4:8]))              # gives the correlation
}

colnames(cand)[9] <- "predictor"
colnames(cand)[10] <- "correlation"

table(cand$predictor)

write.table(cand, file="outliers_correlation.txt", quote = F)

