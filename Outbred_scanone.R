#rQTL script for analyzing an outbred mapping population data

library(qtl)

#Read data into R for csv formated file
data <- read.cross("csv", ".", "Albany_Rqtl.csv", genotypes=NULL)

#Read data into R for a csvr formated file
read.cross("csvr", , "4way_rqtl.csv", genotypes=NULL, alleles=c("A", "B", "C", "D"))

#Jitter markers if necessary
data2<-jittermap(data)

#Imputation of genotype data
dataIM<-calc.genoprob(data2, step=1, map.function="kosambi")

#Check which phenotypes are avaiable for mapping
names(data2$pheno)

#Standard interval mapping using scanone for one trait at a time for Albany
scan=scanone(dataIM, pheno.col="BIOMASS_11.1", model=c("normal"), method=c("hk"))
perm1000=scanone(dataIM, pheno.col="BIOMASS_11.1", model=c("normal"), method=c("hk"), n.perm=1000)
summary(scan, format="tabByCol", pvalues=TRUE, perm=perm1000, ci.function="lodint")
summary(perm1000)
plot(scan)
abline(h=4.33)

#Standard interval mapping using scanone for one trait at a time for new 4way
scan=scanone(dataIM2, pheno.col="Days_Emerge_to_Flower", model=c("normal"), method=c("hk"))
perm1000=scanone(dataIM2, pheno.col="Days_Emerge_to_Flower", model=c("normal"), method=c("hk"), n.perm=1000)
summary(perm1000)
plot(scan)
abline(h=4.47)

#Standard interval mapping using scanone for multiple traits at once 
scan=scanone(dataIM, pheno.col=3:144, model=c("normal"), method=c("hk"))
perm1000=scanone(dataIM, pheno.col=3:144, model=c("normal"), method=c("hk"), n.perm=1000)
tabby<-summary(scan, format="tabByCol", pvalues=TRUE, perm=perm1000, ci.function="lodint")
write.table(tabby, file = "Tab_by_Morphology_2011_means_all.txt", quote=FALSE, row.name=FALSE)

#Estimating allele effects
simdata <- sim.geno(data, step=2, n.draws=128, err=0.001)
qtl <- makeqtl(simdata, chr=18, pos=84.1)
out.fq <- fitqtl(simdata, qtl=qtl, get.ests=TRUE, pheno.col="HRV_MASS_g_12.1", formula=y ~ Q1)

#Visualizing allele effects 
effectplot(data2, pheno.col="Days_Emerge_to_Flower", mname1="10@84")
