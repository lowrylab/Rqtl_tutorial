------------------------------------------------------------------------
Author: David Lowry
Email: davidbryantlowry@gmail.com

This is an example of the workflow for QTL analysis in R/qtl for 
phenotypes measured in the FIL2xHAL2-11 F2 mapping population
------------------------------------------------------------------------


#Loading rQTL library
library(qtl)

#Read in data
data <- read.cross("csv", ".", "Hallii_Rqtl.csv")

#Seperating markers that are in the same location
data<-jittermap(data)

#Building covariate for planting cohort. Only used covariates in final analysis for flowering time.
block1<-as.vector(data$pheno$Cohort1)
block2<-as.vector(data$pheno$Cohort2)
block3<-as.vector(data$pheno$Cohort3)
block4<-as.vector(data$pheno$Cohort4)
x<-cbind(block1, block2, block3, block4)

#Imputation of genotype data
dataIM<-calc.genoprob(data, step=1, map.function="kosambi")

#Check which phenotypes are avaiable for mapping
names(data$pheno)

#Standard interval mapping using scanone 
out2c<-scanone(dataIM, pheno.col=2, addcov=x, model=c("normal"), method=c("hk"))
out2cperm=scanone(dataIM, pheno.col=2, addcov=x, model=c("normal"), method=c("hk"), n.perm=1000)
plot(out2cperm)
summary(out2cperm)
summary(out2c, perms=out2cperm, format="tabByCol", pvalues=TRUE, ci.function="lodint")

#Load modules to allow parallelization across computer cores
require(snow)
require(rlecuyer)

#Loop for scantwo analyses to calaculate penalties for stepwise QTL
#Loop used for multiple phenotypes because scantwo permutations
#are very computationally intensive. This is an example
#where there are 11 phenotypes total

i<-1
while(i<=12){lod2delta_dry <- scantwo(dataIM, pheno.col=i, addcov=x, model=c("normal"), 
method=c("hk")); print(i); print(lod2delta_dry); i<-i+1}

while(i<=12){maxlod2delta <- scantwo(dataIM, pheno.col=i, addcov=x, model=c("normal"),
method=c("hk"), batchsize=100, n.cluster=6, n.perm=1000);pen <-calc.penalties(maxlod2delta); 
print(i); print(pen); i<-i+1}
 
#Stepwise QTL using penalties from scantwo permutations for example of totally additive QTLs 
traitstepwise_11 <- stepwiseqtl(dataIM, pheno.col=11, max.qtl=6, method="hk", penalties=c(3.616570, 5.901286, 3.427231), keeplodprofile=TRUE, keeptrace=TRUE)
plotLodProfile(traitstepwise_11)
plotModel(traitstepwise_11)
summary(traitstepwise_11)

#Calculation of lod scores, variance explained, additive effects, and dominance deviation from best fit stepwise model
morphdataIM_sim <- sim.geno(data, step=1, n.draws=128, err=0.001)
qtl_11 <- makeqtl(morphdataIM_sim, chr=c(1,2,3,7,9), pos=c(83,152,14,26,19))
qtl_11
out.fq_11 <- fitqtl(morphdataIM_sim, qtl=qtl_11, get.ests=TRUE, pheno.col=11, formula=y ~ Q1 + Q2 + Q3 + Q4 + Q5)
out.fq_11

#Stepwise QTL using penalties from scantwo permutations for example of epistatic QTLs. 
traitstepwise_12 <- stepwiseqtl(dataIM, pheno.col=12, max.qtl=6, method="hk", penalties=c(3.631443, 6.495896, 4.187875), keeplodprofile=TRUE, keeptrace=TRUE)

plotLodProfile(traitstepwise_12)
plotModel(traitstepwise_12)
summary(traitstepwise_12)

#Calculation of lod scores, variance explained, additive effects, and dominance deviation from best fit stepwise model
morphdataIM_sim <- sim.geno(data, step=1, n.draws=128, err=0.001)
qtl_12 <- makeqtl(morphdataIM_sim, chr=c(2,7,7), pos=c(65,29,83.1))
qtl_12
out.fq_12 <- fitqtl(morphdataIM_sim, qtl=qtl_12, get.ests=TRUE, pheno.col=2, formula=y ~ Q1 + Q2 + Q3 + Q1:Q3)
out.fq_12

#Calculate 1.5-LOD drop confidence interval for QTLs
lodint(traitstepwise_11, drop=1.5, qtl.index=1, lodcolumn=1)

#Visualizing the effects for one locus
effectplot(dataIM, pheno.col=11, mname1="1@83")

#Visualizing the effects for epistatic interactions for two loci
effectplot(dataIM, pheno.col=12, mname1="2@65", mname2="7@83.1")

