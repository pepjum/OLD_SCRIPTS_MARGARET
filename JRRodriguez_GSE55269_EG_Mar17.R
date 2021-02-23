setwd("/home/bioinformatica/datos/03_Analysis/jrrodriguez/06_GSE55269_Mar17")

# save.image("gse55269.RData")
# load("gse55269.RData")

source("/home/bioinformatica/datos/01_Rscripts/A_Funciones/funcionesVikv2.R")
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

## 8 muestras: silenciamiento y sobre-expresión de SNAI2
gse55269 <- getGEO(filename="GSE55269_family.soft.gz")
charSamples55269 <- list(length(GSMList(gse55269)))
charDescription55269 <- list(length(GSMList(gse55269)))
charTitle55269 <- list(length(GSMList(gse55269)))
for (i in 1:length(GSMList(gse55269))) {                                
	charSamples55269[i] <- list(Meta(GSMList(gse55269)[[i]])$characteristics_ch1)
}
for (i in 1:length(GSMList(gse55269))) {                                
	charDescription55269[i] <- list(Meta(GSMList(gse55269)[[i]])$geo_accession)
}
for (i in 1:length(GSMList(gse55269))) {                                
	charTitle55269[i] <- list(Meta(GSMList(gse55269)[[i]])$title)
}
dataClinic55269 <- data.frame(GEO=unlist(charDescription55269),Sample=unlist(charTitle55269), Cell_type= sapply(strsplit(sapply(charSamples55269,FUN=function(x) x[1]),": "), FUN=function(x) x[2]))
rownames(dataClinic55269) <- dataClinic55269[,1]

# Normalizo pacientes con info sobre KRAS y controles
gse55269_data <- ReadAffy(celfile.path="CEL/")
eset_rma_gse55269 = rma(gse55269_data)
mat_rma_gse55269 = exprs(eset_rma_gse55269)
colnames(mat_rma_gse55269) <- gsub("_HG-U133_Plus_2_.CEL.gz","",sapply(colnames(mat_rma_gse55269),FUN=function(x) unlist(strsplit(x,".Sen."))[2]))
geneNameAll <- lookUp(rownames(mat_rma_gse55269), "hgu133plus2", "SYMBOL");
geneDescAll <- lookUp(rownames(mat_rma_gse55269), "hgu133plus2", "GENENAME");
geneChr <- lookUp(rownames(mat_rma_gse55269), "hgu133plus2", "MAP");
gse55269 = data.frame(Probeset=rownames(mat_rma_gse55269), Name=paste(geneNameAll), Description=paste(geneDescAll), Locus=paste(geneChr), mat_rma_gse55269);
write.table(gse55269, file = "GSE55269_DataRaw.txt", quote = FALSE, row.names = FALSE, sep="\t")

qc.data = qc(gse55269_data)
deg<-AffyRNAdeg(gse55269_data, log.it=T)
acol = sample(brewer.pal(8, "Dark2"), ncol(eset_rma_gse55269), replace = (8 < ncol(eset_rma_gse55269)))
clust.euclid.average <- hclust(dist(t(exprs(gse55269_data))),method="average")
clust.euclid.average2 <- hclust(dist(t(mat_rma_gse55269)),method="average")

pdf(file="GSE55269_QC.pdf", colormode = "rgb", width=20,height=20)
plot(qc.data)
plotAffyRNAdeg(deg, cols = acol, lwd = 2)
hist(gse55269_data)
boxplot(gse55269_data,names=colnames(gse55269_data), cex.axis=0.7, las=2)
plclust(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(as.data.frame(exprs(eset_rma_gse55269)), cex.axis=0.7, las=2)
hist(exprs(eset_rma_gse55269), main="", xlab="Log2(Intendidad)", ylab="Genes")
plclust(clust.euclid.average2, main="Hierarchical clustering of samples (Norm)",  hang=-1)
dev.off()

#estudio variabilidad
type <- factor(c(rep("sh_C",2), rep("sh_SNAI2",2),rep("oe_C",2),rep("oe_SNAI2",2)))
tmp <- detOutliers(mat_rma_gse55269, type, "GSE55269_DetOut.pdf", 16, 16)


## FILTRAR
# 4 grupos

log_intensity_threshold <- 4
numsamples <- 2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.SHC <- genefilter(mat_rma_gse55269[,which(type == "sh_C")], ffun1)
whichFilter.SH <- genefilter(mat_rma_gse55269[,which(type == "sh_SNAI2")], ffun1)
whichFilter.OEC <- genefilter(mat_rma_gse55269[,which(type == "oe_C")], ffun1)
whichFilter.OE <- genefilter(mat_rma_gse55269[,which(type == "oe_SNAI2")], ffun1)
whichFilter <- (whichFilter.SHC | whichFilter.SH | whichFilter.OEC | whichFilter.OE )
sum(whichFilter)

mat_rma_filter <- mat_rma_gse55269[whichFilter,]
mat_rma_filterAnnot <- data.frame(Probeset=rownames(mat_rma_filter), Name=paste(geneNameAll[rownames(mat_rma_filter)]), Description=paste(geneDescAll[rownames(mat_rma_filter)]), Locus=paste(geneChr[rownames(mat_rma_filter)]), mat_rma_filter)
rownames(mat_rma_filterAnnot) <- paste(mat_rma_filterAnnot[,1])
write.table(mat_rma_filterAnnot, file = "DataFilter.txt",quote = FALSE, row.names = FALSE, sep="\t")



#------------------------------------------------------------------------------------------------------------
## LIMMA (4 grupos)

design <- model.matrix(~0+type);
rownames(design) <- colnames(mat_rma_filter)
colnames(design) <- levels(type)
cont.matrix <- makeContrasts(SHvsSHC = sh_SNAI2 - sh_C,
	OEvsOEC = oe_SNAI2 - oe_C,
	levels=design)
fit <- lmFit(as.matrix(mat_rma_filterAnnot[,c(5:ncol(mat_rma_filterAnnot))]), design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

SHvsSHC <- topTable(fit2, coef = "SHvsSHC", n=nrow(fit2), adjust="fdr")
resSHvsSHC <- data.frame(Probeset=rownames(SHvsSHC), Name=paste(geneNameAll[rownames(SHvsSHC)]), Description=paste(geneDescAll[rownames(SHvsSHC)]), Locus=paste(geneChr[rownames(SHvsSHC)]), SHvsSHC)
write.table(resSHvsSHC, file = "SHvsSHC.txt", quote = FALSE, row.names = FALSE, sep="\t")
OEvsOEC <- topTable(fit2, coef = "OEvsOEC", n=nrow(fit2), adjust="fdr")
resOEvsOEC <- data.frame(Probeset=rownames(OEvsOEC), Name=paste(geneNameAll[rownames(OEvsOEC)]), Description=paste(geneDescAll[rownames(OEvsOEC)]), Locus=paste(geneChr[rownames(OEvsOEC)]), OEvsOEC)
write.table(resOEvsOEC, file = "OEvsOEC.txt", quote = FALSE, row.names = FALSE, sep="\t")

pdf(file = "Limma_GSE55269_Mar17.pdf", width = 10, height = 10)
graphContrast(resSHvsSHC, "SH vs SHC (B>0)", 0, 1, 2)
graphContrast(resOEvsOEC, "OE vs OEC (B>0)", 0, 1, 2)
compare3List(paste(resSHvsSHC[resSHvsSHC$B>3,1]), paste(resOEvsOEC[resOEvsOEC$B>3,1]), c(intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC>0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC<0,1])), intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC<0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC>0,1]))) ,"shSNAI2_B3","oeSNAI2_B3", "sh.OE.B3_coherent","Venn Diagram")
compare3List(paste(resSHvsSHC[resSHvsSHC$B>0,1]), paste(resOEvsOEC[resOEvsOEC$B>0,1]), c(intersect(paste(resSHvsSHC[resSHvsSHC$B>0 & resSHvsSHC$logFC>0,1]), paste(resOEvsOEC[resOEvsOEC$B>0 & resOEvsOEC$logFC<0,1])), intersect(paste(resSHvsSHC[resSHvsSHC$B>0 & resSHvsSHC$logFC<0,1]), paste(resOEvsOEC[resOEvsOEC$B>0 & resOEvsOEC$logFC>0,1]))) ,"shSNAI2_B0","oeSNAI2_B0", "sh.OE.B0_coherent","Venn Diagram")
dev.off()

# SNAI2
pdf("SNAI2_Boxplot_Mar17.pdf", width = 10, height = 8)
{ 
	d <- data.frame(as.numeric(mat_rma_filterAnnot[mat_rma_filterAnnot[,2]=="SNAI2",5:ncol(mat_rma_filterAnnot)]), type)
	colnames(d) <- c("expression", "lab")
	boxplot(expression ~ lab, data = d, outline = FALSE, xlab = " ", ylab = "Expression", main = paste("SNAI2 - OE_log2FC=",round(resOEvsOEC[which(resOEvsOEC[,2]=="SNAI2"),"logFC"],3)," (FDR=", round(resOEvsOEC[which(resOEvsOEC[,2]=="SNAI2"),"adj.P.Val"],5), ")", " SH_log2FC=", round(resSHvsSHC[which(resSHvsSHC[,2]=="SNAI2"),"logFC"],3)," (FDR=", round(resSHvsSHC[which(resSHvsSHC[,2]=="SNAI2"),"adj.P.Val"],5), ")",sep=""))
	beeswarm(expression ~ lab, data = d, pch = 21, col = c("blue","orange"), bg = "#00000050", add = TRUE, method="swarm", corral = "random")
}
dev.off()

selcoh_B3 <- mat_rma_filterAnnot[c(intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC>0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC<0,1])), intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC<0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC>0,1]))),]
write.table(cbind(resSHvsSHC[c(intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC>0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC<0,1])), intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC<0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC>0,1]))),], resOEvsOEC[c(intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC>0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC<0,1])), intersect(paste(resSHvsSHC[resSHvsSHC$B>3 & resSHvsSHC$logFC<0,1]), paste(resOEvsOEC[resOEvsOEC$B>3 & resOEvsOEC$logFC>0,1]))),c(5,8:10)]), file = "selcoh_B3.txt", quote = FALSE, row.names = FALSE, sep="\t")

clustData <- as.matrix(selcoh_B3[, c(5:ncol(selcoh_B3))])
rownames(clustData) <- paste(selcoh_B3$Probeset, selcoh_B3$Name)
annotation = data.frame("SampleType" = type)
rownames(annotation) = colnames(clustData)
colType = c("blue", "darkgreen", "darkblue", "darkorange")
names(colType) = c(paste(unique(annotation$SampleType)))#,paste(unique(annotation$Time)))
ann_colors = list("SampleType" = colType)#[1:2], "Time"=colType[3:5])
pheatmap(clustData, annotation = annotation, annotation_colors = ann_colors, border_color = NA, cellwidth = 8, cellheight = 8, fontsize = 8, color = greenred(50), scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cutree_cols = NA, cutree_rows = NA, filename = "Cluster_coherentDEG_B3_Mar17.pdf")

## GO
load("../../vsegura/GOHumanENS_Oct13.RData")
universeRef <- unique(paste(annotEnsembl[paste(annotEnsembl[,2]) %in% paste(unique(paste(mat_rma_filterAnnot$Name))), 1]))

# Coherent DEG (B>3)
ENS_cohB3 <- unique(paste(annotEnsembl[paste(annotEnsembl[,2]) %in% paste(selcoh_B3$Name),1]))
GO_cohB3 <- goGSEnrich(intersect(ENS_cohB3, universeRef), universeRef, gsGOHumanAll, annotEnsembl, 0.01, 2)
GO_cohB3_Filter <- GO_cohB3[which(paste(GO_cohB3[,"GOID"]) %in% names(which(table(GO_cohB3[,"GOID"])>2))),]
write.table(GO_cohB3_Filter, file='GO_cohB3_Filter.txt', quote = FALSE, row.names = TRUE, sep="\t");
GO_cohB3_Filter$logp <- (-1)*log10(as.numeric(paste(GO_cohB3_Filter$pvalue)))
sort1 <- GO_cohB3_Filter[GO_cohB3_Filter[,1] == "BP",]
sort1 <- sort1[order(sort1[, "logp"], decreasing = FALSE), ]
sort2 <- GO_cohB3_Filter[GO_cohB3_Filter[,1] == "MF",]
sort2 <- sort2[order(sort2[, "logp"], decreasing = FALSE), ]
sort3 <- GO_cohB3_Filter[GO_cohB3_Filter[,1] == "CC",]
sort3 <- sort3[order(sort3[, "logp"], decreasing = FALSE), ]
gg1 <- ggplot(data=unique(GO_cohB3_Filter[GO_cohB3_Filter[,1] == "BP",c("GOTerm","logp")]), aes(x = GOTerm, y = logp)) + geom_bar(stat="identity") + labs(x = "", y = "-log(pvalue)") + theme(title=element_text("Biological Process")) + coord_flip() + xlim(unique(paste(sort1[,3]))) 
gg2 <- ggplot(unique(GO_cohB3_Filter[GO_cohB3_Filter[,1] == "MF",c("GOTerm","logp")]), aes(x = GOTerm, y = logp)) + geom_bar(stat="identity") + labs(x = "", y = "-log(pvalue)") + theme(title=element_text("Molecular Function")) + coord_flip() + xlim(unique(paste(sort2[,3]))) 
gg3 <- ggplot(unique(GO_cohB3_Filter[GO_cohB3_Filter[,1] == "CC",c("GOTerm","logp")]), aes(x = GOTerm, y = logp)) + geom_bar(stat="identity") + labs(x = "", y = "-log(pvalue)") + theme( title=element_text("Cellular Component")) + coord_flip() + xlim(unique(paste(sort3[,3]))) 
pdf(file = "GOAnalysis_coherentDEG_B3_Mar17.pdf", height = 35, width = 15, colormodel = "rgb")
gg1
gg2
gg3
dev.off()

