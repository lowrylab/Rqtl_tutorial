------------------------------------------------------------------------
Author: David Lowry
Email: davidbryantlowry@gmail.com

This is an example of the workflow for QTL analysis in R/qtl for 
phenotypes measured in the FIL2xHAL2-11 F2 mapping population
------------------------------------------------------------------------


#Loading rQTL library
library(qtl)

#Read in data
data <- read.cross("csv", ".", "Rqtl_hallii.csv")

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
 
#Stepwise QTL using penalties from scantwo permutations
traitstepwise_2 <- stepwiseqtl(dataIM, pheno.col=2, max.qtl=6, method="hk", 
penalties=c(3.964035, 6.576424, 4.336440), keeplodprofile=TRUE, keeptrace=TRUE)

plotLodProfile(traitstepwise_2)
plotModel(traitstepwise_2)
summary(traitstepwise_2)

#Calculation of lod scores, variance explained, additive effects, and dominance deviation from best fit stepwise model
morphdataIM_sim <- sim.geno(data, step=1, n.draws=128, err=0.001)
qtl_2 <- makeqtl(morphdataIM_sim, chr=c(2,3,5,7), pos=c(48,83,65,58))
qtl_2
out.fq_2 <- fitqtl(morphdataIM_sim, qtl=qtl_2, get.ests=TRUE, pheno.col=2, formula=y ~ Q1 + Q2 + Q3 + Q4 + Q1:Q4)
out.fq_2

#Calculate 1.5-LOD drop confidence interval for QTLs
lodint(traitstepwise_2, drop=1.5, qtl.index=1, lodcolumn=1)

#Visualizing the effects for epistatic interactions
effectplot(morph_data, pheno.col=12, mname1="2@65", mname2="7@83.1")
