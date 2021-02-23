setwd("/home/bioinformatica/datos/03_Analysis/eguruce/16_aaindex_Mar16")

# save.image("aaindex_mar16.RData")
# load("aaindex_mar16.RData")


###########################################################################
# 1. Pruebo los paquetes para estudiar las propiedades de los péptidos
# 2. Selecciono péptidos proteotípicos de NextProt con observaciones en GPMDB y péptidos proteotípicos missing sin observaciones en GPMDB
# 3. Calculo la media de cada propiedad para los péptidos seleccionados
# 4. Comparo con t-test MissingObs vs NoMissingObs y MissingNObs vs NoMissing Obs
# 5. Preparo matriz para hacer resampling y sacar un p-value robusto para cada propiedad
# 6. Selección de propiedades y péptidos para el clasificador en "aaindex_classifier2_EG_Abr16.R" (No sirve)
# 6b. Estudio los péptidos más o menos observados en GPMDB sin tener en cuenta si son missing "aaindex_gpmdb_EG_may"



# Uso el paquete seqinr de CRAN que tiene la versión 9.1 de aaindex
# install.packages("/home/calculus/03_Analysis/eguruce/16_aaindex_Mar16/seqinr_3.1-3.tar.gz", repos=NULL, type="source")
# install.packages("/home/calculus/03_Analysis/eguruce/16_aaindex_Mar16/Interpol_1.3.1.tar.gz", repos=NULL, type="source")
# install.packages("/home/calculus/03_Analysis/eguruce/16_aaindex_Mar16/protr_1.1-1.tar.gz", repos=NULL, type="source")

library(seqinr)
data(aaindex) # seqinr::aaindex (544 properties - aaindex version 9.1)
library(protr) # protr::AAindex (544 properties) podría intentar usarse con la función de Interpol)
library(Interpol) # normaliza las properties pero utiliza AAindex y no aaindex (Interpol::AAindex (533 properties)) [AAdescriptor(data, descriptor = 151, normalize = 0)], interpola la secuencia de una proteina a un tamaño determinado para el clasificador [Interpol(data, dims, method = "linear")]
library(doBy)

AAindex_protr <- AAindex
AAindex <- AAindex[,-c(2:6)]

# seqinr
which(sapply(aaindex, function(x) length(grep("Size", x$D)) != 0))
a()
aaa()
a("Ala")
aaa("A")

ex <- c("MALRGKRKAEDF", c2s(rev(s2c("MALRGKRKAEDF"))))
names(ex) <- c("ex1","ex2")

stats <- AAstat(s2c("MALRGKRKAEDF"))
computePI(s2c("MALRGKRKAEDF"))
getLength("MALRGKRKAEDF") # 12
pmw(s2c("MALRGKRKAEDF"))
getName()
getTrans()
compoAA <- table(factor(s2c(ex[1]), levels = LETTERS))
nTermR <- which(LETTERS == s2c(ex[1])[1])
cTermR <- which(LETTERS == s2c(ex[1])[nchar(ex[1])])

dotPlot(s2c("MALRGKRKAE"),s2c("MALRGKRKAE"), main="Internal repeats")
dotPlot(s2c("MALRGKRKAE"),rev(s2c("MALRGKRKAE")), main="Inversion")

# Interpol
aaindex[[2]]$I[aaa(s2c("MALRGKRKAEDF"))]
AAdescriptor("MALRGKRKAEDF",2,0)
hydrophobicity_norm <- AAdescriptor("MALRGKRKAEDF",2,1)
hydrophobicity_norm2 <- AAdescriptor("MALRGKRKAEDF",2,2) # Entre 0 y 1

hydrophobicity_norm_l10 <- Interpol(hydrophobicity_norm2, 10, method = "linear")
plot(unlist(hydrophobicity_norm2), type="l", col="darkgreen", ylim=c(-3,3), ylab="descriptor value", xlab="sequence position", lwd=2)
lines(seq(1,12,12/10), as.vector(hydrophobicity_norm_l10),col="red",lwd=2)
axis(3, at=seq(1,12,12/10), labels=1:10)

# protR
protcheck(ex)

extractAPAAC(paste0(ex[1],ex[2],ex[1])) # con proteínas de más de 30 aa o cambiar lambda
# AA composition
extractAAC(ex[1]) # Amino acid composition
extractDC(ex[2]) # Dipeptide composition
extractTC(ex[1]) # Tripeptide composition

# Protein similarity
twoSeqSim(ex[1], ex[2])
parSeqSim(list(ex[1],ex[2]))
twoGOSim()
parGOSim()

# Autocorrelation
extractMoreauBroto(paste0(ex[1],ex[2],ex[1]))
extractMoran(ex[1])
extractGeary()
# CTD
extractCTDC(ex[1]) # Composition
extractCTDT(ex[1]) # Transition
extractCTDD(ex[1]) # Distribution
# Conjoint triad descriptors
extractCTriad(ex[1])
# Quasi-sequence-order descriptors
extractSOCN(paste0(ex[1],ex[2],ex[1])) # Sequence-order-coupling number
extractQSO(paste0(ex[1],ex[2],ex[1])) # Quasi-sequence-order descriptors
# Pseudo-amino acid composition
extractPAAC(paste0(ex[1],ex[2],ex[1]))
extractAPAAC(paste0(ex[1],ex[2],ex[1]))
# Profile-based descriptors
extractPSSM(ex[1])
extractPSSMAcc()
extractPSSMFeature()


#----------------------------------------------------------------------------------------------------------------------------------------------------
# PEPTIDES (GPMDB and NextProt)


# PEPTIDOS proteotípicos de NEXTPROT versión 01/2016 (script "\\valis\CALCULUS\01_Rscripts\05_AnalisisProteomicaMS/Analisis_MissingProteins_PRIDE_NextProt20150901_AG_Ene16.R" en "\\Valis\calculus\03_Analysis\agarin\NextProt_20160111\")
library(gdata)
NextProtXChr <- data.frame()
chrList <- c(1:22, "X", "Y", "MT")
for (i in 1:length(chrList)) {
	cat(i, "\n")
	filetmp <- paste("/home/calculus/03_Analysis/agarin/NextProt_20160111/nextprot_chromosome_", chrList[i], ".txt", sep = "")
	con <- file(filetmp, "r", blocking = FALSE)
	geneListv2 <- readLines(con)
	geneListv2.lst <- strsplit(geneListv2[21:(length(geneListv2)-7)], "\t")
	geneListv2.tmp <- t(as.data.frame(lapply(geneListv2.lst, FUN = function(x) x[[1]])))
	geneListv2.PE <- trim(t(as.data.frame(lapply(geneListv2.lst, FUN = function(x) x[[3]]))))
	geneListv2.tmp2 <- strsplit(paste(geneListv2.tmp[,1])," ")
	geneListv2.tmp2 <- lapply(geneListv2.tmp2, FUN = function(x) x[x != ""])

	if (i == 1) {
		nextProtXChr <- data.frame(chrList[i], sapply(geneListv2.tmp2, FUN = function(x) x[[length(x)]]), sapply(geneListv2.tmp2, FUN = function(x) x[[1]]), geneListv2.PE)
	} else {
		nextProtXChr <- rbind(nextProtXChr, data.frame(chrList[i], sapply(geneListv2.tmp2, FUN = function(x) x[[length(x)]]), sapply(geneListv2.tmp2, FUN = function(x) x[[1]]), geneListv2.PE))
	}
}
rownames(nextProtXChr) <- NULL
colnames(nextProtXChr) <- c("Chr", "NextprotID", "GeneName", "ProteinEvidence")
tmp <- paste(nextProtXChr$ProteinEvidence)
tmp[tmp == "protein level"] = "PE1"
tmp[tmp == "transcript level"] = "PE2"
tmp[tmp == "homology"] = "PE3"
tmp[tmp == "predicted"] = "PE4"
tmp[tmp == "uncertain"] = "PE5"
nextProtXChr$PE <- tmp
nextProtXChr$Missing <- (nextProtXChr$PE %in% c("PE2", "PE3", "PE4"))*1
save(nextProtXChr, file = "nextProtXChr_20160111.RData")


fileName <- "nextProtDB20160111_simple.txt"
con <- file(fileName, open="r")
MissingPepAll <- readLines(con)
close(con)
protsAll <- MissingPepAll[grep("Protein Name", MissingPepAll)]
protAnnotAll <- t(data.frame(strsplit(unlist(lapply(strsplit(protsAll, ": "), FUN = function(x) x[2])), "\\|")))
colnames(protAnnotAll) <- NULL
rownames(protAnnotAll) <- NULL
npeptidesAll <- MissingPepAll[grep("Number", MissingPepAll)]
NpeptidesAll <- as.numeric(paste(unlist(lapply(strsplit(npeptidesAll, " = "), FUN = function(x) x[2]))))
tmp <- MissingPepAll[-grep("Protein Name", MissingPepAll)]
tmp <- tmp[tmp != ""]
tmp <- tmp[-grep("Number", tmp)]
tmp <- tmp[-grep("No", tmp)]
tmp <- tmp[-c(1,2)]
tmp2 <- strsplit(tmp, " ")
tmp2 <- lapply(tmp2, FUN = function(x) x[x != ""])
peptidesAll.df <- data.frame("PepNo" = unlist(lapply(tmp2, FUN = function(x) x[1])), "Range" = unlist(lapply(tmp2, FUN = function(x) x[2])), "IsotopicMass" = unlist(lapply(tmp2, FUN = function(x) x[3])), "AverageMass" = unlist(lapply(tmp2, FUN = function(x) x[4])), "Peptide" = unlist(lapply(tmp2, FUN = function(x) x[5])))  # 7031728 proteotypic pep - 2958550 unique pep
npepind <- c()
for (i in 1:length(NpeptidesAll))
	npepind <- c(npepind, rep(i, NpeptidesAll[i]))
peptidesXProtAll.df <- data.frame(protAnnotAll[npepind,], peptidesAll.df)
colnames(peptidesXProtAll.df)[1] <- "NextprotID"
peptidesXProtAll.df$Range <- paste("'", peptidesXProtAll.df$Range, sep = "")
peptidesXProtAll.df$Prot <- sapply(paste(peptidesXProtAll.df[,1]), FUN=function(x) unlist(strsplit(x,"-"))[1])
peptidesXProtAll.df <- merge(peptidesXProtAll.df, nextProtXChr[,c("NextprotID","Missing")],by.x=7, by.y=1,all.x=T, all.y=F)
pepXmissing <- summaryBy(Missing~Peptide, data=peptidesXProtAll.df[,c("Peptide","Missing")], fun=min)
pepXmissing[pepXmissing==Inf] <- NA

# Descarga de GPMDB el 04/04/2016 ("http://peptides.thegpm.org/~/peptides_by_species/")
# Péptidos de NextProt vistos en GPMDB
pepObs_gpmdb <- read.table("Homo_sapiens.tsv", sep="\t", header=T) # 1195063 peptides
pepObs_gpmdb_proteotypic <- pepObs_gpmdb[pepObs_gpmdb[,1] %in% paste(unique(peptidesAll.df[,"Peptide"])),]
pepObs_gpmdb_proteotypicXmis <- merge(pepObs_gpmdb_proteotypic, pepXmissing, by.x=1, by.y=1, all.x=T, all.y=F)
# Filtro los NextProt_ID de cromosomas raros para los que tengo péptidos pero no información de si es missing
pepObs_gpmdb_proteotypicXmis_f <- pepObs_gpmdb_proteotypicXmis[!is.na(pepObs_gpmdb_proteotypicXmis$Missing.min),]
pepObs_gpmdb_proteotypicXmis_f2 <- pepObs_gpmdb_proteotypicXmis_f[sapply(paste(pepObs_gpmdb_proteotypicXmis_f[,1]), FUN=function(x) protcheck(x)),]
save("pepObs_gpmdb_proteotypicXmis_f2",file="processedPepInput.rda")
# Péptidos de NextProt no vistos en GPMDB y que sean missing (casos negativos para el clasificador)
pepNObs_proteotypic <- unique(pepXmissing[(!(paste(pepXmissing[,1]) %in% paste(pepObs_gpmdb[,1]))) & (pepXmissing[,2]==1),])
pepNObs_proteotypic_f <- pepNObs_proteotypic[!is.na(pepNObs_proteotypic$Missing.min),]
pepNObs_proteotypic_f2 <- pepNObs_proteotypic_f[sapply(paste(pepNObs_proteotypic_f[,1]), FUN=function(x) protcheck(x)),]
save("pepNObs_proteotypic_f2",file="processedPepNObsInput.rda")

#--------------------------------------------------------------------------------------------------------------------------------------------------
# PROPERTIES per Peptide
# "For each peptide and a given property, the constituent amino acid numerical values were averaged to produce a single value. Missing values were ignored. The average (rather than median or sum) was chosen because it is sensitive to outliers and normalizes for peptide length. It was assumed that the average physicochemical property across each peptide was sufficient to capture relevant information about peptide response." Fusaro et al. 2009

# prop_annot <- sapply(aaindex, FUN=function(x) x$D)
# size_mean <- mean(aaindex[[names(which(sapply(aaindex, function(x) length(grep("Size", x$D)) != 0)))]]$I[aaa(s2c("MALRGKRKAEDF"))])
# prop_mean <- c(lapply(aaindex, FUN=function(x) mean(x$I[aaa(s2c("MALRGKRKAEDF"))], na.rm=TRUE)), LENGTH=getLength("MALRGKRKAEDF"), PMW=pmw(s2c("MALRGKRKAEDF")), PI=computePI(s2c("MALRGKRKAEDF")), stats$Prop)

# OBS
pep_length <- sapply(paste(pepObs_gpmdb_proteotypicXmis_f2[,1]), FUN=function(x) getLength(x))
pep_PMW <-  sapply(paste(pepObs_gpmdb_proteotypicXmis_f2[,1]), FUN=function(x) pmw(s2c(x)))
pep_PI <- sapply(paste(pepObs_gpmdb_proteotypicXmis_f2[,1]), FUN=function(x) computePI(s2c(x)))
# N-OBS
pep_length_nobs <- sapply(paste(pepNObs_proteotypic_f2[,1]), FUN=function(x) getLength(x))


# En el cluster
# pep_aaStat <- sapply(paste(pepObs_gpmdb_proteotypicXmis_f2[,1]), FUN=function(x) AAstat(s2c(x), plot=FALSE)$Prop)
load("1_pep_aaStat.rda")
prop_mean <- list()
# prop_mean <- sapply(paste(pepObs_gpmdb_proteotypicXmis_f2[,1]), FUN=function(x) lapply(aaindex, FUN=function(y) mean(y$I[aaa(s2c(x))], na.rm=TRUE)), LENGTH=getLength(x), PMW=pmw(s2c(x)), PI=computePI(s2c(x)), AAstat(s2c(x), plot=FALSE)$Prop))

# Ejecutado en paralelo en el cluster
# prop_mean <- sapply(paste(pepObs_gpmdb_proteotypicXmis_f2[1:10,1]), FUN=function(x) lapply(aaindex, FUN=function(y) mean(y$I[aaa(s2c(x))], na.rm=TRUE)))
load("aaindex_avg1.RDA")
prop_mean_c <- prop_mean
load("aaindex_avg2.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg3.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg4.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg5.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg6.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg7.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg8.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg9.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg10.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg11.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg12.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean)
load("aaindex_avg13.RDA")
prop_mean_c <- cbind(prop_mean_c, prop_mean[,1:18])
rm(prop_mean)

all_prop <- rbind(prop_mean_c, pep_aaStat, Length=pep_length, PMW=pep_PMW, PI=pep_PI)
save("all_prop", file="all_prop.rda")
for(i in 1:10)
{
	all_prop_c <- all_prop[c((i*50-49):(i*50)),]
	save(all_prop_c, file=paste(i,"_all_prop.rda",sep=""))
}
all_prop_c <- all_prop[501:nrow(all_prop),]
save(all_prop_c, file="11_all_prop.rda")

# NObs
load("1_pepNObs_aaindex_AVG.rda")
propNObs_mean_c <- prop_mean
load("50001_pepNObs_aaindex_AVG.rda")
propNObs_mean_c <- cbind(propNObs_mean_c, prop_mean)
load("100001_pepNObs_aaindex_AVG.rda")
propNObs_mean_c <- cbind(propNObs_mean_c, prop_mean)
load("150001_pepNObs_aaindex_AVG.rda")
propNObs_mean_c <- cbind(propNObs_mean_c, prop_mean)
rm(prop_mean)

load("1_pepNObs_aaStat.rda")
pepNObs_aaStat_c <- pep_aaStat
load("50001_pepNObs_aaStat.rda")
pepNObs_aaStat_c <- cbind(pepNObs_aaStat_c, pep_aaStat)
load("100001_pepNObs_aaStat.rda")
pepNObs_aaStat_c <- cbind(pepNObs_aaStat_c, pep_aaStat)
load("150001_pepNObs_aaStat.rda")
pepNObs_aaStat_c <- cbind(pepNObs_aaStat_c, pep_aaStat)
load("200001_pepNObs_aaStat.rda")
pepNObs_aaStat_c <- cbind(pepNObs_aaStat_c, pep_aaStat)
rm(pep_aaStat)

#pep_PMW_nobs <-  sapply(paste(pepNObs_proteotypic_f2[,1]), FUN=function(x) pmw(s2c(x)))
load("pwm/1_pepNObs_PMW.rda")
pep_PMW_nobs_c <- pep_PMW
load("pwm/50001_pepNObs_PMW.rda")
pep_PMW_nobs_c <- c(pep_PMW_nobs_c, pep_PMW)
load("pwm/100001_pepNObs_PMW.rda")
pep_PMW_nobs_c <- c(pep_PMW_nobs_c, pep_PMW)
load("pwm/150001_pepNObs_PMW.rda")
pep_PMW_nobs_c <- c(pep_PMW_nobs_c, pep_PMW)
rm(pep_PMW)

# pep_PI_nobs <- sapply(paste(pepNObs_proteotypic_f2[,1]), FUN=function(x) computePI(s2c(x)))
load("pi/1_pepNObs_PI.RDA")
pep_PI_nobs_c <- pep_PI
load("pi/2_pepNObs_PI.RDA")
pep_PI_nobs_c <- c(pep_PI_nobs_c, pep_PI)
load("pi/3_pepNObs_PI.RDA")
pep_PI_nobs_c <- c(pep_PI_nobs_c, pep_PI)
load("pi/4_pepNObs_PI.RDA")
pep_PI_nobs_c <- c(pep_PI_nobs_c, pep_PI)
rm(pep_PI)

all_propNObs <- rbind(propNObs_mean_c, pepNObs_aaStat_c, Length=pep_length_nobs, PMW=pep_PMW_nobs_c, PI=pep_PI_nobs_c)
save("all_propNObs", file="all_propNObs.rda")
for(i in 1:10)
{
	all_propNObs_c <- all_propNObs[c((i*50-49):(i*50)),]
	save(all_propNObs_c, file=paste(i,"_all_propNObs.rda",sep=""))
}
all_propNObs_c <- all_propNObs[501:nrow(all_propNObs),]
save(all_propNObs_c, file="11_all_propNObs.rda")

#--------------------------------------------------------------------------------------------------------------------
## T-test:
# Missing_gpmdb vs NoMissing_gpmdb
L_MissVsNoMiss_t <- t.test(pep_length[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1],pep_length[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0], na.rm=TRUE)
PMW_MissVsNoMiss_t <- t.test(pep_PMW[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1],pep_PMW[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0], na.rm=TRUE)
PI_MissVsNoMiss_t <- t.test(pep_PI[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1],pep_PI[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0], na.rm=TRUE)
aaStat_MissVsNoMiss_t <- apply(pep_aaStat, 1, FUN=function(x) t.test(unlist(x)[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1],unlist(x)[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0], na.rm=TRUE))
Prop_MissVsNoMiss_t <- apply(prop_mean_c, 1, FUN=function(x) t.test(unlist(x)[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1],unlist(x)[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0], na.rm=TRUE))
# Missing_NoGpmdb vs NoMissing_gpmdb
L_MissNObsVsNoMiss_t <- t.test(pep_length_nobs,unlist(all_prop["Length",pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0]), na.rm=TRUE)
PMW_MissNObsVsNoMiss_t <- t.test(pep_PMW_nobs_c,unlist(all_prop["PMW",pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0]), na.rm=TRUE)
PI_MissNObsVsNoMiss_t <- t.test(pep_PI_nobs_c,unlist(all_prop["PI",pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0]), na.rm=TRUE)
aaStat_MissNObsVsNoMiss_t <- sapply(1:nrow(pepNObs_aaStat_c), FUN=function(x) t.test(unlist(pepNObs_aaStat_c[as.numeric(paste(x)),]),unlist(all_prop[as.numeric(paste(x))+544,pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0]), na.rm=TRUE))
Prop_MissNObsVsNoMiss_t <- sapply(1:nrow(propNObs_mean_c),FUN=function(x) t.test(unlist(propNObs_mean_c[as.numeric(paste(x)),]),unlist(all_prop[as.numeric(paste(x)),pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0]), na.rm=TRUE))

# RDA con los t-test
save(list=c("L_MissVsNoMiss_t","PMW_MissVsNoMiss_t","PI_MissVsNoMiss_t","aaStat_MissVsNoMiss_t","Prop_MissVsNoMiss_t","L_MissNObsVsNoMiss_t","PMW_MissNObsVsNoMiss_t","PI_MissNObsVsNoMiss_t","aaStat_MissNObsVsNoMiss_t","Prop_MissNObsVsNoMiss_t"), file="global_ttest_res.rda")

#--------------------------------------------------------------------------------------------------------------
## "Rarefaction" para comprobar que el p-value es robusto
# 500 iteraciones con el 80 % de los genes missing (rownames->prop , colnames=iteration)
sample_500 <- cbind(sapply(c(1:500),FUN=function(x) sample(which(pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1),round(sum(pepObs_gpmdb_proteotypicXmis_f2$Missing.min)*0.8))), sapply(c(1:500),FUN=function(x) sample(which(pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0),round(sum(pepObs_gpmdb_proteotypicXmis_f2$Missing.min)*0.8))))
save("sample_500", file="sample_500.rda")

# En el cluster
# rarefaction_pvalue <- sapply(c(1:500), FUN=function(x) apply(all_prop, 1, FUN=function(y) t.test(unlist(y)[sample_500[,as.numeric(paste(x))]],unlist(y)[sample_500[,as.numeric(paste(x))+500]], na.rm=TRUE)$p.value))
load(".rda")
rarefaction_pvalue_c <- rarefaction_pvalue

# Control -: 500 iteraciones con el 80 % de los genes no-missing vs 80% de genes no-missing(rownames->prop , colnames=iteration)
sample_500_Control <- cbind(sapply(c(1:500),FUN=function(x) sample(which(pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0),round(sum(pepObs_gpmdb_proteotypicXmis_f2$Missing.min)*0.8))), sapply(c(1:500),FUN=function(x) sample(which(pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0),round(sum(pepObs_gpmdb_proteotypicXmis_f2$Missing.min)*0.8))))
save("sample_500_Control", file="sample_500_Control.rda")

# En el cluster
# rarefaction_pvalue_controls <- sapply(c(1:500), FUN=function(x) apply(all_prop[1:3,], 1, FUN=function(y) t.test(unlist(y)[sample_500_Control[,as.numeric(paste(x))]],unlist(y)[sample_500_Control[,as.numeric(paste(x))+500]], na.rm=TRUE)$p.value))

# 500 iteraciones con el 80 % de los genes missing y no observados (rownames->prop , colnames=iteration)
sampleNObs_500 <- cbind(sapply(c(1:500),FUN=function(x) sample(1:nrow(pepNObs_proteotypic_f2),round(nrow(pepNObs_proteotypic_f2)*0.8))), sapply(c(1:500),FUN=function(x) sample(which(pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0),round(nrow(pepNObs_proteotypic_f2)*0.8))))
save("sampleNObs_500", file="sampleNObs_500.rda")
