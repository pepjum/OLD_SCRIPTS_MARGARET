setwd("/home/nostromo/data/pepe/60_Arrays_TESIS_UPEFINDER_PRRG3/GSE4290/")

source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
library(affy)
library(simpleaffy)
library(RColorBrewer)
library(limma)
library(doBy)
library(GEOquery)
library(hgu133plus2.db)
library(pheatmap)
library(gplots)
library(annotate)
library(beeswarm)
library(ggplot2)



## 23 samples from epilepsy patients are used as non-tumor samples. 157 tumor samples include 26 astrocytomas, 50 oligodendrogliomas and 81 glioblastomas. =180
gse4290 <- getGEO(filename="GSE4290_family.soft.gz")
charSamples4290 <- list(length(GSMList(gse4290)))
charDescription4290 <- list(length(GSMList(gse4290)))
charTitle4290 <- list(length(GSMList(gse4290)))
for (i in 1:length(GSMList(gse4290))) {                                
	charSamples4290[i] <- list(Meta(GSMList(gse4290)[[i]])$characteristics_ch1)
}
for (i in 1:length(GSMList(gse4290))) {                                
	charDescription4290[i] <- list(Meta(GSMList(gse4290)[[i]])$geo_accession)
}
for (i in 1:length(GSMList(gse4290))) {                                
	charTitle4290[i] <- list(Meta(GSMList(gse4290)[[i]])$title)
}
dataClinic4290 <- data.frame(GEO=unlist(charDescription4290),Sample=unlist(charTitle4290), Cell_type= sapply(strsplit(sapply(charSamples4290,FUN=function(x) x[1]),": "), FUN=function(x) x[2]))
rownames(dataClinic4290) <- dataClinic4290[,1]

dataClinic4290$Cell_type<-gsub(" ","_",dataClinic4290$Cell_type)
dataClinic4290$Cell_type[which(is.na(dataClinic4290$Cell_type))]<-"tumor"


##unzip GZ files
files_to_unzip<-list.files("../GSE4290_RAW/", pattern=".gz")
files_to_unzip<-paste0(getwd(),"_RAW/",files_to_unzip)

for(file in files_to_unzip){
	y<-paste0("gunzip -d ", file)
	system(y)
}

tmp2<-dataClinic4290[c("GEO","Cell_type")]

gse4290_data <- ReadAffy(celfile.path="../GSE4290_RAW/")
eset_rma_gse4290 = rma(gse4290_data)
mat_rma_gse4290 = exprs(eset_rma_gse4290)
#colnames(mat_rma_gse4290) <- gsub("_HG-U133_Plus_2_.CEL.gz","",sapply(colnames(mat_rma_gse4290),FUN=function(x) unlist(strsplit(x,".Sen."))[2]))
geneNameAll <- lookUp(rownames(mat_rma_gse4290), "hgu133plus2", "SYMBOL");
geneDescAll <- lookUp(rownames(mat_rma_gse4290), "hgu133plus2", "GENENAME");
geneChr <- lookUp(rownames(mat_rma_gse4290), "hgu133plus2", "MAP");
gse4290 = data.frame(Probeset=rownames(mat_rma_gse4290), Name=paste(geneNameAll), Description=paste(geneDescAll), Locus=paste(geneChr), mat_rma_gse4290);
write.table(gse4290, file = "GSE4290_DataRaw.txt", quote = FALSE, row.names = FALSE, sep="\t")

qc.data = qc(gse4290_data)
deg<-AffyRNAdeg(gse4290_data, log.it=T)
acol = sample(brewer.pal(180, "Dark2"), ncol(eset_rma_gse4290), replace = (8 < ncol(eset_rma_gse4290)))
clust.euclid.average <- hclust(dist(t(exprs(gse4290_data))),method="average")
clust.euclid.average2 <- hclust(dist(t(mat_rma_gse4290)),method="average")


pdf(file="GSE4290_QC.pdf", colormode = "rgb", width=20,height=20)
plot(qc.data)
plotAffyRNAdeg(deg, cols = acol, lwd = 2)
hist(gse4290_data)
boxplot(gse4290_data,names=colnames(gse4290_data), cex.axis=0.7, las=2)
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(as.data.frame(exprs(eset_rma_gse4290)), cex.axis=0.7, las=2)
hist(exprs(eset_rma_gse4290), main="", xlab="Log2(Intendidad)", ylab="Genes")
plot(clust.euclid.average2, main="Hierarchical clustering of samples (Norm)",  hang=-1)
dev.off()


#outliers

#type<-rownames(gse4290)
mat_rma_gse4290_2<-mat_rma_gse4290
colnames(mat_rma_gse4290_2)<-paste0(lapply(strsplit(paste(colnames(mat_rma_gse4290_2)),".CEL"),"[",1),"_",tmp2$Cell_type)
colnames(mat_rma_gse4290_2)<-tmp2$Cell_type

type<-factor(colnames(mat_rma_gse4290_2))
type<-str_replace(type,",","")
type<-str_replace(type,"-","")


tmp <- detOutliers(mat_rma_gse4290_2, type, "GSE4290_DetOut_2.pdf", 16, 16)


## FILTRAR
# 7 grupos

log_intensity_threshold <- 5

numsamples <- 3   #ajustar para que sea el 50% de muestras de cada tipo
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)
whichFilter.a2 <- genefilter(mat_rma_gse4290_2[,which(type == "astrocytoma_grade_2")], ffun1)

numsamples <- 9   #ajustar para que sea el 50% de muestras de cada tipo
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)
whichFilter.a3 <- genefilter(mat_rma_gse4290_2[,which(type == "astrocytoma_grade_3")], ffun1)

numsamples <- 38   #ajustar para que sea el 50% de muestras de cada tipo
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)
whichFilter.g4 <- genefilter(mat_rma_gse4290_2[,which(type == "glioblastoma_grade_4")], ffun1)

numsamples <- 11   #ajustar para que sea el 50% de muestras de cada tipo
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)
whichFilter.c <- genefilter(mat_rma_gse4290_2[,which(type == "nontumor")], ffun1)

numsamples <- 2   #ajustar para que sea el 50% de muestras de cada tipo
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)
whichFilter.t <- genefilter(mat_rma_gse4290_2[,which(type == "tumor")], ffun1)

numsamples <- 19   #ajustar para que sea el 50% de muestras de cada tipo
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)
whichFilter.o2 <- genefilter(mat_rma_gse4290_2[,which(type == "oligodendroglioma_grade_2")], ffun1)

numsamples <- 6   #ajustar para que sea el 50% de muestras de cada tipo
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)
whichFilter.o3 <- genefilter(mat_rma_gse4290_2[,which(type == "oligodendroglioma_grade_3")], ffun1)


whichFilter <- (whichFilter.a2 | whichFilter.a3 | whichFilter.g4 | whichFilter.c | whichFilter.t | whichFilter.o2 | whichFilter.o3 )
sum(whichFilter)  #31440	

mat_rma_filter <- mat_rma_gse4290_2[whichFilter,]
mat_rma_filterAnnot <- data.frame(Probeset=rownames(mat_rma_filter), Name=paste(geneNameAll[rownames(mat_rma_filter)]), Description=paste(geneDescAll[rownames(mat_rma_filter)]), Locus=paste(geneChr[rownames(mat_rma_filter)]), mat_rma_filter)
rownames(mat_rma_filterAnnot) <- paste(mat_rma_filterAnnot[,1])
write.table(mat_rma_filterAnnot, file = "DataFilter.txt",quote = FALSE, row.names = FALSE, sep="\t")


#------------------------------------------------------------------------------------------------------------
## LIMMA (7 grupos)
library(stringr)
colnames(mat_rma_filter)<-str_replace(colnames(mat_rma_filter),",","")
colnames(mat_rma_filter)<-str_replace(colnames(mat_rma_filter),"-","")

type<-str_replace(type,",","")
type<-str_replace(type,"-","")

type<-factor(type)
design <- model.matrix(~0+type);
rownames(design) <- colnames(mat_rma_filter)
colnames(design) <- levels(type)
cont.matrix <- makeContrasts(a2vsC = astrocytoma_grade_2 - nontumor,
							 a3vsC = astrocytoma_grade_3 - nontumor,
							 tvsC = tumor - nontumor,
							 gvsC = glioblastoma_grade_4 - nontumor,
							 o2vsC = oligodendroglioma_grade_2 - nontumor,
							 o3vsC = oligodendroglioma_grade_3 - nontumor,
	levels=design)


#ajustar al modelo lineal
fit <- lmFit(as.matrix(mat_rma_filterAnnot[,c(5:ncol(mat_rma_filterAnnot))]), design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2) #empiricalBayesforDifferencial expression

a2vsC <- topTable(fit2, coef = "a2vsC", n=nrow(fit2), adjust="fdr")
resa2vsC <- data.frame(Probeset=rownames(a2vsC), Name=paste(geneNameAll[rownames(a2vsC)]), Description=paste(geneDescAll[rownames(a2vsC)]), Locus=paste(geneChr[rownames(a2vsC)]), a2vsC)
write.table(resa2vsC, file = "resa2vsC.txt", quote = FALSE, row.names = FALSE, sep="\t")

a3vsC <- topTable(fit2, coef = "a3vsC", n=nrow(fit2), adjust="fdr")
resa3vsC <- data.frame(Probeset=rownames(a3vsC), Name=paste(geneNameAll[rownames(a3vsC)]), Description=paste(geneDescAll[rownames(a3vsC)]), Locus=paste(geneChr[rownames(a3vsC)]), a3vsC)
write.table(resa3vsC, file = "resa3vsC", quote = FALSE, row.names = FALSE, sep="\t")

tvsC <- topTable(fit2, coef = "tvsC", n=nrow(fit2), adjust="fdr")
restvsC <- data.frame(Probeset=rownames(tvsC), Name=paste(geneNameAll[rownames(tvsC)]), Description=paste(geneDescAll[rownames(tvsC)]), Locus=paste(geneChr[rownames(tvsC)]), tvsC)
write.table(restvsC, file = "restvsC", quote = FALSE, row.names = FALSE, sep="\t")

gvsC <- topTable(fit2, coef = "gvsC", n=nrow(fit2), adjust="fdr")
resgvsC <- data.frame(Probeset=rownames(gvsC), Name=paste(geneNameAll[rownames(gvsC)]), Description=paste(geneDescAll[rownames(gvsC)]), Locus=paste(geneChr[rownames(gvsC)]), gvsC)
write.table(resgvsC, file = "resgvsC", quote = FALSE, row.names = FALSE, sep="\t")

o2vsC <- topTable(fit2, coef = "o2vsC", n=nrow(fit2), adjust="fdr")
reso2vsC <- data.frame(Probeset=rownames(o2vsC), Name=paste(geneNameAll[rownames(o2vsC)]), Description=paste(geneDescAll[rownames(o2vsC)]), Locus=paste(geneChr[rownames(o2vsC)]), o2vsC)
write.table(reso2vsC, file = "reso2vsC", quote = FALSE, row.names = FALSE, sep="\t")

o3vsC <- topTable(fit2, coef = "o3vsC", n=nrow(fit2), adjust="fdr")
reso3vsC <- data.frame(Probeset=rownames(o3vsC), Name=paste(geneNameAll[rownames(o3vsC)]), Description=paste(geneDescAll[rownames(o3vsC)]), Locus=paste(geneChr[rownames(o3vsC)]), o3vsC)
write.table(reso3vsC, file = "reso3vsC", quote = FALSE, row.names = FALSE, sep="\t")


pdf(file = "Limma_GSE4290_NOV20.pdf", width = 10, height = 10)
graphContrast(resa2vsC, "A2 vs C (B>0)", 0, 1, 2)
graphContrast(resa3vsC, "A3 vs C (B>0)", 0, 1, 2)
graphContrast(restvsC, "T vs C (B>0)", 0, 1, 2)
graphContrast(resgvsC, "Glio vs C (B>0)", 0, 1, 2)
graphContrast(reso2vsC, "O2 vs C (B>0)", 0, 1, 2)
graphContrast(reso3vsC, "O3 vs C (B>0)", 0, 1, 2)
#compare3List(paste(resSHvsSHC[resSHvsSHC$B>3,1]), paste(resOEvsOEC[resOEvsOEC$B>3,1]), c(intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC>0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC<0,1])), intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC<0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC>0,1]))) ,"shSNAI2_B3","oeSNAI2_B3", "sh.OE.B3_coherent","Venn Diagram")
#compare3List(paste(resSHvsSHC[resSHvsSHC$B>0,1]), paste(resOEvsOEC[resOEvsOEC$B>0,1]), c(intersect(paste(resSHvsSHC[resSHvsSHC$B>0 & resSHvsSHC$logFC>0,1]), paste(resOEvsOEC[resOEvsOEC$B>0 & resOEvsOEC$logFC<0,1])), intersect(paste(resSHvsSHC[resSHvsSHC$B>0 & resSHvsSHC$logFC<0,1]), paste(resOEvsOEC[resOEvsOEC$B>0 & resOEvsOEC$logFC>0,1]))) ,"shSNAI2_B0","oeSNAI2_B0", "sh.OE.B0_coherent","Venn Diagram")
dev.off()


# PRRG3
#separar por muestras   y por canales si hay mas de uno para ese gen
d <- data.frame(mat_rma_filterAnnot[which(mat_rma_filterAnnot[,2]=="PRRG3"),][5:ncol(mat_rma_filterAnnot)])
names_astro_grade3<-c("astrocytoma._grade_3",paste0("astrocytoma._grade_3.", seq(1:100)))
names_astro_grade2<-c("astrocytoma._grade_2",paste0("astrocytoma._grade_2.", seq(1:100)))
names_nontumor<-c("non.tumor",paste0("non.tumor.", seq(1:100)))
names_glio4<-c("glioblastoma._grade_4",paste0("glioblastoma._grade_4.", seq(1:100)))
names_oligo2<-c("oligodendroglioma._grade_2",paste0("oligodendroglioma._grade_2.", seq(1:100)))
names_oligo3<-c("oligodendroglioma._grade_3",paste0("oligodendroglioma._grade_3.", seq(1:100)))


nontumor_df<- d[,colnames(d) %in% c(names_nontumor,"Probeset")]
nontumor_df_ch1<-nontumor_df[c("220433_at"),]
nontumor_df_ch2<-nontumor_df[c("229118_at"),]

glio_df<-d[,colnames(d) %in% c(names_glio4,"Probeset")]
glio_df_ch1<-glio_df[c("220433_at"),]
glio_df_ch2<-glio_df[c("229118_at"),]


as2_df<-d[,colnames(d) %in% c(names_astro_grade2,"Probeset")]
as2_df_ch1<-as2_df[c("220433_at"),]
as2_df_ch2<-as2_df[c("229118_at"),]

as3_df<-d[,colnames(d) %in% c(names_astro_grade3,"Probeset")]
as3_df_ch1<-as3_df[c("220433_at"),]
as3_df_ch2<-as3_df[c("229118_at"),]

ol2_df<-d[,colnames(d) %in% c(names_oligo2,"Probeset")]
ol2_df_ch1<-ol2_df[c("220433_at"),]
ol2_df_ch2<-ol2_df[c("229118_at"),]

ol3_df<-d[,colnames(d) %in% c(names_oligo3,"Probeset")]
ol3_df_ch1<-ol3_df[c("220433_at"),]
ol3_df_ch2<-ol3_df[c("229118_at"),]


list_nont_ch1<-c()
for(i in 1:ncol(nontumor_df_ch1)){
	k<-nontumor_df_ch1[,i]
	list_nont_ch1<-c(list_nont_ch1,k)
}

df_nont_ch1<-data.frame("values"=list_nont_ch1,"sample"=rep("nontumor",length(list_nont_ch1)))

list_nont_ch2<-c()
for(i in 1:ncol(nontumor_df_ch2)){
	k<-nontumor_df_ch2[,i]
	list_nont_ch2<-c(list_nont_ch2,k)
}

df_nont_ch2<-data.frame("values"=list_nont_ch2,"sample"=rep("nontumor",length(list_nont_ch2)))


list_glio_ch1<-c()
for(i in 1:ncol(glio_df_ch1)){
	k<-glio_df_ch1[,i]
	list_glio_ch1<-c(list_glio_ch1,k)
}

df_glio_ch1<-data.frame("values"=list_glio_ch1,"sample"=rep("gliog4",length(list_glio_ch1)))

list_glio_ch2<-c()
for(i in 1:ncol(glio_df_ch2)){
	k<-glio_df_ch2[,i]
	list_glio_ch2<-c(list_glio_ch2,k)
}

df_glio_ch2<-data.frame("values"=list_glio_ch2,"sample"=rep("gliog4",length(list_glio_ch2)))


list_as2_ch1<-c()
for(i in 1:ncol(as2_df_ch1)){
	k<-as2_df_ch1[,i]
	list_as2_ch1<-c(list_as2_ch1,k)
}

df_as2_ch1<-data.frame("values"=list_as2_ch1,"sample"=rep("as2",length(list_as2_ch1)))

list_as2_ch2<-c()
for(i in 1:ncol(as2_df_ch2)){
	k<-as2_df_ch2[,i]
	list_as2_ch2<-c(list_as2_ch2,k)
}

df_as2_ch2<-data.frame("values"=list_as2_ch2,"sample"=rep("as2",length(list_as2_ch2)))


list_as3_ch1<-c()
for(i in 1:ncol(as3_df_ch1)){
	k<-as3_df_ch1[,i]
	list_as3_ch1<-c(list_as3_ch1,k)
}

df_as3_ch1<-data.frame("values"=list_as3_ch1,"sample"=rep("as3",length(list_as3_ch1)))

list_as3_ch2<-c()
for(i in 1:ncol(as3_df_ch2)){
	k<-as3_df_ch2[,i]
	list_as3_ch2<-c(list_as3_ch2,k)
}

df_as3_ch2<-data.frame("values"=list_as3_ch2,"sample"=rep("as3",length(list_as3_ch2)))


list_ol2_ch1<-c()
for(i in 1:ncol(ol2_df_ch1)){
	k<-ol2_df_ch1[,i]
	list_ol2_ch1<-c(list_ol2_ch1,k)
}

df_ol2_ch1<-data.frame("values"=list_ol2_ch1,"sample"=rep("ol2",length(list_ol2_ch1)))

list_ol2_ch2<-c()
for(i in 1:ncol(ol2_df_ch2)){
	k<-ol2_df_ch2[,i]
	list_ol2_ch2<-c(list_ol2_ch2,k)
}

df_ol2_ch2<-data.frame("values"=list_ol2_ch2,"sample"=rep("ol2",length(list_ol2_ch2)))

list_ol3_ch1<-c()
for(i in 1:ncol(ol3_df_ch1)){
	k<-ol3_df_ch1[,i]
	list_ol3_ch1<-c(list_ol3_ch1,k)
}

df_ol3_ch1<-data.frame("values"=list_ol3_ch1,"sample"=rep("ol3",length(list_ol3_ch1)))

list_ol3_ch2<-c()
for(i in 1:ncol(ol3_df_ch2)){
	k<-ol3_df_ch2[,i]
	list_ol3_ch2<-c(list_ol3_ch2,k)
}

df_ol3_ch2<-data.frame("values"=list_ol3_ch2,"sample"=rep("ol3",length(list_ol3_ch2)))

plotter1<-rbind(df_nont_ch1, df_glio_ch1,df_as2_ch1, df_as3_ch1, df_ol2_ch1, df_ol3_ch1 )
plotter2<-rbind(df_nont_ch2, df_glio_ch2,df_as2_ch2, df_as3_ch2, df_ol2_ch2, df_ol3_ch2 )



pdf("PRRG3_Boxplot_Nov20_ch1.pdf", width = 10, height = 8)
boxplot(values ~ sample, data=plotter1, outline=FALSE, xlab="220433_at", ylab="Expression", main=c('','',paste0("pval glio4vsC= ", signif(resgvsC[which(resgvsC$Name=="PRRG3" & resgvsC$Probeset=="220433_at"), "P.Value"],2)),paste0("pval as2vsC= ",signif(resa2vsC[which(resa2vsC$Name=="PRRG3" & resa2vsC$Probeset=="220433_at"), "P.Value"],2)), paste0("pval as3vsC= ",signif(resa3vsC[which(resa3vsC$Name=="PRRG3" & resa3vsC$Probeset=="220433_at"), "P.Value"],2)), paste0("pval ol2vsC= ",signif(reso2vsC[which(reso2vsC$Name=="PRRG3" & reso2vsC$Probeset=="220433_at"), "P.Value"],2)),paste0("pval ol3vsC= ",signif(reso3vsC[which(reso3vsC$Name=="PRRG3" & reso3vsC$Probeset=="220433_at"), "P.Value"],2))),cex.main=1)
dev.off()

pdf("PRRG3_Boxplot_Nov20_ch1_ggplo2.pdf", width = 10, height = 8)
pl<-ggplot(plotter1, aes(x=sample, y=values))
pl + geom_boxplot(aes(fill=sample)) + theme_bw() + ggtitle("PRRG3 en canal 220433_at")
dev.off()



pdf("PRRG3_Boxplot_Nov20_ch2.pdf", width = 10, height = 8)
boxplot(values ~ sample, data=plotter2, outline=FALSE, xlab="229118_at", ylab="Expression", main=c('','',paste0("pval glio4vsC= ", signif(resgvsC[which(resgvsC$Name=="PRRG3" & resgvsC$Probeset=="229118_at"), "P.Value"],2)),paste0("pval as2vsC= ",signif(resa2vsC[which(resa2vsC$Name=="PRRG3" & resa2vsC$Probeset=="229118_at"), "P.Value"],2)), paste0("pval as3vsC= ",signif(resa3vsC[which(resa3vsC$Name=="PRRG3" & resa3vsC$Probeset=="229118_at"), "P.Value"],2)), paste0("pval ol2vsC= ",signif(reso2vsC[which(reso2vsC$Name=="PRRG3" & reso2vsC$Probeset=="229118_at"), "P.Value"],2)),paste0("pval ol3vsC= ",signif(reso3vsC[which(reso3vsC$Name=="PRRG3" & reso3vsC$Probeset=="229118_at"), "P.Value"],2))),cex.main=1)
dev.off()

pdf("PRRG3_Boxplot_Nov20_ch2_ggplo2.pdf", width = 10, height = 8)
pl<-ggplot(plotter2, aes(x=sample, y=values))
pl + geom_boxplot(aes(fill=sample)) + theme_bw() + ggtitle("PRRG3 en canal 229118_at")
dev.off()



###### crear df_con las columnas de pvalues para mi gen

reso2vsC_sel<-reso2vsC[which(reso2vsC$Name=="PRRG3"),c("Probeset","Name","P.Value","t","B")]
names(reso2vsC_sel)[3:5]<-paste0("O2vsC_",names(reso2vsC_sel)[3:5])
reso3vsC_sel<-reso3vsC[which(reso3vsC$Name=="PRRG3"),c("Probeset","Name","P.Value","t","B")]
names(reso3vsC_sel)[3:5]<-paste0("O3vsC_",names(reso3vsC_sel)[3:5])
resgvsC_sel<-resgvsC[which(resgvsC$Name=="PRRG3"),c("Probeset","Name","P.Value","t","B")]
names(resgvsC_sel)[3:5]<-paste0("GlvsC_",names(resgvsC_sel)[3:5])
resa2vsC_sel<-resa2vsC[which(resa2vsC$Name=="PRRG3"),c("Probeset","Name","P.Value","t","B")]
names(resa2vsC_sel)[3:5]<-paste0("a2vsC_",names(resa2vsC_sel)[3:5])
resa3vsC_sel<-resa3vsC[which(resa3vsC$Name=="PRRG3"),c("Probeset","Name","P.Value","t","B")]
names(resa3vsC_sel)[3:5]<-paste0("a3vsC_",names(resa3vsC_sel)[3:5])

all_PRRG3<-cbind(reso2vsC_sel,reso3vsC_sel[,3:5],resgvsC_sel[,3:5],resa2vsC_sel[,3:5],resa3vsC_sel[,3:5])   # esto esta bien para ver los pvalues y las B en cada muestra



#### seleccionar un grupo de genes de cada muestra para ver funciones anotadas en GO

##PRRG3 B en reso2vsC = 16.417314

reso2vsC_sel_ALL<-reso2vsC[which(reso2vsC$B> 16),]  ### 742 Genes

resgvsC_sel_ALL<-resgvsC[which(resgvsC$B> 31.2),] ### 994 genes

## EXTRAER ENSG ASOCIADO A CADA Genename
load("~/data/03_Analysis/vsegura/GOHumanENS_Oct13.RData")
universeRef <- unique(paste(annotEnsembl[paste(annotEnsembl[,2]) %in% paste(unique(paste(mat_rma_filterAnnot$Name))), 1]))

ENSG_resgvC_sel_ALL <- unique(paste(annotEnsembl[paste(annotEnsembl[,2]) %in% paste(unique(paste(resgvsC_sel_ALL$Name))), 1]))

GO_cohB30 <- goGSEnrich(intersect(ENSG_resgvC_sel_ALL, universeRef), universeRef, gsGOHumanAll, annotEnsembl, 0.01, 2)
GO_cohB30_Filter <- GO_cohB30[which(paste(GO_cohB30[,"GOID"]) %in% names(which(table(GO_cohB30[,"GOID"])>2))),]

write.table(GO_cohB30_Filter, file='GO_cohB30_Filter.txt', quote = FALSE, row.names = TRUE, sep="\t");
GO_cohB30_Filter$logp <- (-1)*log10(as.numeric(paste(GO_cohB30_Filter$pvalue)))

sort1 <- GO_cohB30_Filter[GO_cohB30_Filter[,1] == "BP",]
sort1 <- sort1[order(sort1[, "logp"], decreasing = FALSE), ]
sort1_sel<-sort1[which(sort1$logp >6),]

sort2 <- GO_cohB30_Filter[GO_cohB30_Filter[,1] == "MF",]
sort2 <- sort2[order(sort2[, "logp"], decreasing = FALSE), ]
sort2_sel<-sort2[which(sort2$logp >4),]


sort3 <- GO_cohB30_Filter[GO_cohB30_Filter[,1] == "CC",]
sort3 <- sort3[order(sort3[, "logp"], decreasing = FALSE), ]
sort3_sel<-sort3[which(sort3$logp >4),]


blue.bold.italic.0.01.text <- element_text(face = "bold.italic", color = "blue", size = 20)

pdf('GOAnalysis_coherent_DEG_B30_NOV20_BP_best.pdf', 40,40)
 ggplot(data=unique(GO_cohB30_Filter[GO_cohB30_Filter[,1] == "BP",c("GOTerm","logp")]), aes(x = GOTerm, y = logp)) + geom_bar(stat="identity") + labs(x = "", y = "-log(pvalue)") + theme(axis.text.y = blue.bold.italic.0.01.text) + ggtitle("Biological Process")+ coord_flip() + xlim(unique(paste(sort1_sel[,3]))) 

dev.off()
pdf('GOAnalysis_coherent_DEG_B30_NOV20_MF_best.pdf', 40,40)
ggplot(unique(GO_cohB30_Filter[GO_cohB30_Filter[,1] == "MF",c("GOTerm","logp")]), aes(x = GOTerm, y = logp)) + geom_bar(stat="identity") + labs(x = "", y = "-log(pvalue)") + theme(axis.text.y = blue.bold.italic.0.01.text) + ggtitle("Molecular Function") + coord_flip() + xlim(unique(paste(sort2_sel[,3]))) 
dev.off()

pdf('GOAnalysis_coherent_DEG_B30_NOV20_CC_best.pdf', 40,40)
ggplot(unique(GO_cohB30_Filter[GO_cohB30_Filter[,1] == "CC",c("GOTerm","logp")]), aes(x = GOTerm, y = logp)) + geom_bar(stat="identity") + labs(x = "", y = "-log(pvalue)") + theme(axis.text.y = blue.bold.italic.0.01.text) + ggtitle("Cellular Component") + coord_flip() + xlim(unique(paste(sort3_sel[,3]))) 
dev.off()


######### 10 controles y 12 cancer colorrectal  GSE4107
dir.create("/home/nostromo/data/pepe/60_Arrays_TESIS_UPEFINDER_PRRG3/GSE4107/")
setwd("/home/nostromo/data/pepe/60_Arrays_TESIS_UPEFINDER_PRRG3/GSE4107/")

source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
library(affy)
library(simpleaffy)
library(RColorBrewer)
library(limma)
library(doBy)
library(GEOquery)
library(hgu133plus2.db)
library(pheatmap)
library(gplots)
library(annotate)
library(beeswarm)
library(ggplot2)



gse4107 <- getGEO(filename="GSE4107_family.soft.gz")
charSamples4107 <- list(length(GSMList(gse4107)))
charDescription4107 <- list(length(GSMList(gse4107)))
charTitle4107 <- list(length(GSMList(gse4107)))
for (i in 1:length(GSMList(gse4107))) {                                
	charSamples4107[i] <- list(Meta(GSMList(gse4107)[[i]])$characteristics_ch1)
}
for (i in 1:length(GSMList(gse4107))) {                                
	charDescription4107[i] <- list(Meta(GSMList(gse4107)[[i]])$geo_accession)
}
for (i in 1:length(GSMList(gse4107))) {                                
	charTitle4107[i] <- list(Meta(GSMList(gse4107)[[i]])$title)
}
dataClinic4107 <- data.frame(GEO=unlist(charDescription4107),Sample=unlist(charTitle4107), Cell_type= paste(charSamples4107))
rownames(dataClinic4107) <- dataClinic4107[,1]
dataClinic4107$Sample<-lapply(strsplit(paste(dataClinic4107$Sample)," "),"[",3)

files_to_unzip<-list.files("../GSE4107_RAW/", pattern=".gz")
files_to_unzip<-paste0(getwd(),"_RAW/",files_to_unzip)

for(file in files_to_unzip){
	y<-paste0("gunzip -d ", file)
	system(y)
}





gse4107_data <- ReadAffy(celfile.path="/home/nostromo/data/pepe/60_Arrays_TESIS_UPEFINDER_PRRG3/GSE4107_RAW/")
eset_rma_gse4107 = rma(gse4107_data)
mat_rma_gse4107 = exprs(eset_rma_gse4107)
colnames(mat_rma_gse4107) <- gsub(".CEL","", colnames(mat_rma_gse4107))
colnames(mat_rma_gse4107)<-paste0(colnames(mat_rma_gse4107),"_",dataClinic4107$Sample)
geneNameAll <- lookUp(rownames(mat_rma_gse4107), "hgu133plus2", "SYMBOL");
geneDescAll <- lookUp(rownames(mat_rma_gse4107), "hgu133plus2", "GENENAME");
geneChr <- lookUp(rownames(mat_rma_gse4107), "hgu133plus2", "MAP");
gse4107 = data.frame(Probeset=rownames(mat_rma_gse4107), Name=paste(geneNameAll), Description=paste(geneDescAll), Locus=paste(geneChr), mat_rma_gse4107);
write.table(gse4107, file = "GSE4107_DataRaw.txt", quote = FALSE, row.names = FALSE, sep="\t")

qc.data = qc(gse4107_data)
deg<-AffyRNAdeg(gse4107_data, log.it=T)
acol = sample(brewer.pal(8, "Dark2"), ncol(eset_rma_gse4107), replace = (8 < ncol(eset_rma_gse4107)))
clust.euclid.average <- hclust(dist(t(exprs(gse4107_data))),method="average")
clust.euclid.average2 <- hclust(dist(t(mat_rma_gse4107)),method="average")

pdf(file="GSE4107_QC.pdf", colormode = "rgb", width=20,height=20)
plot(qc.data)
plotAffyRNAdeg(deg, cols = acol, lwd = 2)
hist(gse4107_data)
boxplot(gse4107_data,names=colnames(gse4107_data), cex.axis=0.7, las=2)
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(as.data.frame(exprs(eset_rma_gse4107)), cex.axis=0.7, las=2)
hist(exprs(eset_rma_gse4107), main="", xlab="Log2(Intendidad)", ylab="Genes")
plot(clust.euclid.average2, main="Hierarchical clustering of samples (Norm)",  hang=-1)
dev.off()


#estudio variabilidad
type<-factor(c(rep("ccolo",12), rep("control",10)))

tmp <- detOutliers(mat_rma_gse4107, type, "GSE4107_DetOut.pdf", 16, 16)


log_intensity_threshold <- 4
numsamples <- 4
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.ccolo <- genefilter(mat_rma_gse4107[,which(type == "ccolo")], ffun1)
whichFilter.control <- genefilter(mat_rma_gse4107[,which(type == "control")], ffun1)
whichFilter <- (whichFilter.ccolo | whichFilter.control)
sum(whichFilter) #22682

mat_rma_filter <- mat_rma_gse4107[whichFilter,]
mat_rma_filterAnnot <- data.frame(Probeset=rownames(mat_rma_filter), Name=paste(geneNameAll[rownames(mat_rma_filter)]), Description=paste(geneDescAll[rownames(mat_rma_filter)]), Locus=paste(geneChr[rownames(mat_rma_filter)]), mat_rma_filter)
rownames(mat_rma_filterAnnot) <- paste(mat_rma_filterAnnot[,1])
write.table(mat_rma_filterAnnot, file = "DataFilter.txt",quote = FALSE, row.names = FALSE, sep="\t")


## LIMMA (4 grupos)

design <- model.matrix(~0+type);
rownames(design) <- colnames(mat_rma_filter)
colnames(design) <- levels(type)
cont.matrix <- makeContrasts(CcovsC = ccolo - control,
	levels=design)

fit <- lmFit(as.matrix(mat_rma_filterAnnot[,c(5:ncol(mat_rma_filterAnnot))]), design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

CcovsC <- topTable(fit2, coef = "CcovsC", n=nrow(fit2), adjust="fdr")
resCcovsC <- data.frame(Probeset=rownames(CcovsC), Name=paste(geneNameAll[rownames(CcovsC)]), Description=paste(geneDescAll[rownames(CcovsC)]), Locus=paste(geneChr[rownames(CcovsC)]), CcovsC)
write.table(resCcovsC, file = "CcovsC.txt", quote = FALSE, row.names = FALSE, sep="\t")

##### NOSALE



