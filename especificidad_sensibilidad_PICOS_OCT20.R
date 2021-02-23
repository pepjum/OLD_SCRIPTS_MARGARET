
genecode<-"/home/nostromo/data/00_References/gencode_gtf/gencode.v34.annotation.gtf"

library(ChIPpeakAnno)


parseENCODE <- function(x) {

	tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
	tmp <- tmp[tmp != ""]
	tmp <- tmp[seq(2, 12, 2)]
	return(tmp)
}

### AHORA CARGAMOS GENCODE y SELECCIONULLMOS LAS COLUMNULLS QUE NOS INTERESAN PARA ESTE TIPO DE ANULLLISIS


genecodev25_all <- read.table(genecode, skip = 5, header = FALSE, sep = "\t")

genecodev25 <- genecodev25_all[genecodev25_all$V3 == "gene", ]
genecodev25_tmp <- apply(as.data.frame(genecodev25[,9]), 1, parseENCODE)
genecodev25_tmp <- t(genecodev25_tmp)

colnames(genecodev25_tmp) <- c("gene_id", "gene_type", "gene_status", "gene_NULLme", "level", "havaNULL_gene")
genecodev25_Annot <- data.frame(genecodev25[,1:8], genecodev25_tmp[, c(1, 2,3, 4, 5,6)])

genecodev25_Annot_Gene <- genecodev25_Annot

G25AnnotChIP <- unique(genecodev25_Annot[genecodev25_Annot[, "V3"] == "gene", c("V1", "V4", "V5", "V7", "gene_id")])
colNULLmes(G25AnnotChIP) <- c("chr", "start", "end", "strand", "gene_id")
G25AnnotChIP_RL <- RangedData(IRanges(start = G25AnnotChIP$start, end = G25AnnotChIP$end, NULLmes = paste(G25AnnotChIP$gene_id)), space = paste(G25AnnotChIP$chr), strand = paste(G25AnnotChIP$strand))


selection<-genecodev25_Annot[,c("gene_status","gene_id")]

#### pasar a ENSG la lista de genes de MSIGDB

EGR1_MSIG_DB<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/EGR1_M1532.txt", skip=2)

EGR1_with_ENSG<-merge(EGR1_MSIG_DB,selection, by.x="V1", by.y="gene_status")
write.table(EGR1_with_ENSG, file="/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/EGR1_M1532_with_ENSG.txt", col.NULLmes=F, row.NULLmes=F, quote=F, sep="\t")

SP1_MSIG_DB<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1_M4831.txt", skip=2)

SP1_with_ENSG<-merge(SP1_MSIG_DB,selection, by.x="V1", by.y="gene_status")
write.table(SP1_with_ENSG, file="/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1_M1532_with_ENSG.txt", col.NULLmes=F, row.NULLmes=F, quote=F, sep="\t")

JUND_MSIG_DB<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND_M7937.txt", skip=2)

JUND_with_ENSG<-merge(JUND_MSIG_DB,selection, by.x="V1", by.y="gene_status")
write.table(JUND_with_ENSG, file="/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND_M1532_with_ENSG.txt", col.NULLmes=F, row.NULLmes=F, quote=F, sep="\t")


#### METRICAS EGR1

#ensembl<-60669
#msigdb_EGR1<-270 ENSG

results_EGR1<-data.frame("METHOD"=NA,"TP"=NA,"TN"=NA,"FP"=NA,"FN"=NA,"TPR"=NA,"TNR"=NA,"ACC"=NA,"PPV"=NA,"F1"=NA)

PEAKZCL_EGR1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/EGR1/PEAKZCL/PeakZCL_peaks_zCros_22_clus_1_norm_SES_min_100_max_1500.sorted_final_annot.txt", header=T)

TP<-length(intersect(paste(PEAKZCL_EGR1$feature),paste(EGR1_with_ENSG$gene_id)))
TN<-60669-length(unique(paste(PEAKZCL_EGR1$feature)))
FN<-270-TP
FP<-length(unique(paste(PEAKZCL_EGR1$feature)))-270
TPR<-TP/(TP+FN)
TNR<-TN/(TN+TP)
ACC<-(TP+TN)/(TP+TN+FP+FN)
PPV<-TP/(TP+FP)
F1<-2*(TP)/(2*(TP)+FP+FN)

results_EGR1<-rbind(results_EGR1,c("PEAKZCL", TP, TN, FP, FN,TPR,TNR,ACC,PPV,F1))
results_EGR1<-results_EGR1[-c(1),]
MACS2_EGR1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/EGR1/MACS/MACS2_EGR1_peaks_annot.txt", header=T)

TP<-length(intersect(paste(MACS2_EGR1$feature),paste(EGR1_with_ENSG$gene_id)))
TN<-60669-length(unique(paste(MACS2_EGR1$feature)))
FN<-270-TP
FP<-length(unique(paste(MACS2_EGR1$feature)))-270
TPR<-TP/(TP+FN)
TNR<-TN/(TN+TP)
ACC<-(TP+TN)/(TP+TN+FP+FN)
PPV<-TP/(TP+FP)
F1<-2*(TP)/(2*(TP)+FP+FN)

results_EGR1<-rbind(results_EGR1,c("MACS2", TP, TN, FP, FN,TPR,TNR,ACC,PPV,F1))

RANGER_EGR1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/EGR1/RANGER/RANGER_EGR1.bed_region_annot.txt", header=T)

TP<-length(intersect(paste(RANGER_EGR1$feature),paste(EGR1_with_ENSG$gene_id)))
TN<-60669-length(unique(paste(RANGER_EGR1$feature)))
FN<-270-TP
FP<-length(unique(paste(RANGER_EGR1$feature)))-270
TPR<-TP/(TP+FN)
TNR<-TN/(TN+TP)
ACC<-(TP+TN)/(TP+TN+FP+FN)
PPV<-TP/(TP+FP)
F1<-2*(TP)/(2*(TP)+FP+FN)

results_EGR1<-rbind(results_EGR1,c("RANGER", TP, TN, FP, FN,TPR,TNR,ACC,PPV,F1))

JAMM_EGR1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/EGR1/JAMM/selected_JAMM_EGR1_annot.txt", header=T)

TP<-length(intersect(paste(JAMM_EGR1$feature),paste(EGR1_with_ENSG$gene_id)))
TN<-60669-length(unique(paste(JAMM_EGR1$feature)))
FN<-270-TP
FP<-length(unique(paste(JAMM_EGR1$feature)))-270
TPR<-TP/(TP+FN)
TNR<-TN/(TN+TP)
ACC<-(TP+TN)/(TP+TN+FP+FN)
PPV<-TP/(TP+FP)
F1<-2*(TP)/(2*(TP)+FP+FN)

results_EGR1<-rbind(results_EGR1,c("JAMM", TP, TN, FP, FN,TPR,TNR,ACC,PPV,F1))


SICER_EGR1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/EGR1/SICER/CHIP_sorted_EGR1-W200-G200-islands-summary-FDR0.01_annot.txt", header=T)

TP<-length(intersect(paste(SICER_EGR1$feature),paste(EGR1_with_ENSG$gene_id)))
TN<-60669-length(unique(paste(SICER_EGR1$feature)))
FN<-270-TP
FP<-length(unique(paste(SICER_EGR1$feature)))-270
TPR<-TP/(TP+FN)
TNR<-TN/(TN+TP)
ACC<-(TP+TN)/(TP+TN+FP+FN)
PPV<-TP/(TP+FP)
F1<-2*(TP)/(2*(TP)+FP+FN)

results_EGR1<-rbind(results_EGR1,c("SICER", TP, TN, FP, FN,TPR,TNR,ACC,PPV,F1))


BAYES_EGR1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/EGR1/BAYES/BayesPeak_output_EGR1_annot.txt", header=T)

TP<-length(intersect(paste(BAYES_EGR1$feature),paste(EGR1_with_ENSG$gene_id)))
TN<-60669-length(unique(paste(BAYES_EGR1$feature)))
FN<-270-TP
FP<-length(unique(paste(BAYES_EGR1$feature)))-270
TPR<-TP/(TP+FN)
TNR<-TN/(TN+TP)
ACC<-(TP+TN)/(TP+TN+FP+FN)
PPV<-TP/(TP+FP)
F1<-2*(TP)/(2*(TP)+FP+FN)

results_EGR1<-rbind(results_EGR1,c("BAYES", TP, TN, FP, FN,TPR,TNR,ACC,PPV,F1))

#### no se habla de esto en el paper. Solo calculado para H3k4me3 que tenemos un ground truth bueno. Para los factores de trascripcion solo hicimos phyper y MEME

phyper_ZCL<-phyper(overlap, group2, Total-group2,group1, lower.tail=T)
overlap=TP
group2=total dianas MSIGDB
total-group2<-total genes-total dianas MSIGDB
group1<-numero de picos que da el algoritmo

hiper_ZCL<-phyper(114-1,270,(60669-270),10697,lower.tail = FALSE) ### 3.313525e-21
hiper_MACS<-phyper(172-1, 270,(60669-270), 10766, lower.tail=F) ## 5.331287e-63
hiper_RANGER<-phyper(185-1, 270,(60669-270), 14319, lower.tail=F) ## 3.883395e-55
hiper_JAMM<-phyper(190-1, 270,(60669-270), 15825, lower.tail=F)  #2.494958e-52
hiper_SICER<-phyper(185-1, 270,(60669-270), 16248, lower.tail=F) # 1.768755e-46
hiper_BAYES<-phyper(208-1, 270,(60669-270), 24290, lower.tail=F) # 2.694368e-35

##### hipergeometrica de JUND


MACS2_JUND<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND/MACS/MACS2_JUND_peaks_annot.txt", header=T)
RANGER_JUND<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND/RANGER/RANGER_JUND.bed_region_annot.txt", header=T)
SICER_JUND<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND/SICER/CHIP_sorted_JUND-W200-G200-islands-summary-FDR0.01_annot.txt", header=T)
BAYES_JUND<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND/BAYES/BayesPeak_output_JUND_annot.txt", header=T)
JAMM_JUND<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND/JAMM/selected_JAMM_JUND_annot.txt", header=T)
PEAKZCL_JUND<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND/PEAKZCL/PeakZCL_peaks_zCros_24_clus_1_norm_SES_min_60_max_1000.sorted_final_annot.txt", header=T)


JUND_MSIG_DB<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND_M7937.txt", skip=2)

JUND_with_ENSG<-merge(JUND_MSIG_DB,selection, by.x="V1", by.y="gene_status")
write.table(JUND_with_ENSG, file="/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/JUND_M7937_with_ENSG.txt", col.names=F, row.names=F, quote=F, sep="\t")

TP<-length(intersect(paste(MACS2_JUND$feature),paste(JUND_with_ENSG$gene_id)))
hiper_MACS2<-phyper(170-1, 270,(60669-270), 17899, lower.tail=F) ## 4.647049e-30

TP<-length(intersect(paste(RANGER_JUND$feature),paste(JUND_with_ENSG$gene_id)))
hiper_RANGER<-phyper(190-1, 270,(60669-270), 22701, lower.tail=F) ## 3.864335e-28

TP<-length(intersect(paste(SICER_JUND$feature),paste(JUND_with_ENSG$gene_id)))
hiper_SICER<-phyper(182-1, 270,(60669-270), 23648, lower.tail=F) ## 2.861971e-21

TP<-length(intersect(paste(JAMM_JUND$feature),paste(JUND_with_ENSG$gene_id)))
hiper_JAMM<-phyper(188-1, 270,(60669-270), 22446, lower.tail=F) ## 1.222019e-27

TP<-length(intersect(paste(PEAKZCL_JUND$feature),paste(JUND_with_ENSG$gene_id)))
hiper_PEAKZCL<-phyper(100-1, 270,(60669-270), 12634, lower.tail=F) ## 6.588657e-10



M4831 


MACS2_SP1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1/MACS/MACS2_SP1_peaks_annot.txt", header=T)
RANGER_SP1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1/RANGER/RANGER_SP1.bed_region_annot.txt", header=T)
SICER_SP1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1/SICER/CHIP_sorted_SP1-W200-G200-islands-summary-FDR0.01_annot.txt", header=T)
BAYES_SP1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1/BAYES/BayesPeak_output_SP1_annot.txt", header=T)
JAMM_SP1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1/JAMM/selected_JAMM_SP1_annot.txt", header=T)
PEAKZCL_SP1<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1/PEAKZCL/PeakZCL_peaks_zCros_24_clus_1_norm_SES_min_120_max_1000.sorted_final_annot.txt", header=T)


SP1_MSIG_DB<-read.table("/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1_M4831.txt", skip=2)

SP1_with_ENSG<-merge(SP1_MSIG_DB,selection, by.x="V1", by.y="gene_status")
write.table(SP1_with_ENSG, file="/home/nostromo/data/pepe/RESULTADOS_PEAKZCL_OCT2020/SP1_M4831_with_ENSG.txt", col.names=F, row.names=F, quote=F, sep="\t")

TP<-length(intersect(paste(MACS2_SP1$feature),paste(SP1_with_ENSG$gene_id)))
hiper_MACS2<-phyper(88-1, 239,(60669-239), 5463, lower.tail=F) ## 5.210283e-32

TP<-length(intersect(paste(RANGER_SP1$feature),paste(SP1_with_ENSG$gene_id)))
hiper_RANGER<-phyper(80-1, 239,(60669-239), 6064, lower.tail=F) ## 3.964605e-23

TP<-length(intersect(paste(SICER_SP1$feature),paste(SP1_with_ENSG$gene_id)))
hiper_SICER<-phyper(139-1, 270,(60669-270), 16273, lower.tail=F) ## 6.418252e-18

TP<-length(intersect(paste(JAMM_SP1$feature),paste(SP1_with_ENSG$gene_id)))
hiper_JAMM<-phyper(126-1, 270,(60669-270), 14178, lower.tail=F) ## 4.083874e-17

TP<-length(intersect(paste(PEAKZCL_SP1$feature),paste(SP1_with_ENSG$gene_id)))
hiper_PEAKZCL<-phyper(41-1, 270,(60669-270), 4135, lower.tail=F) ## 1.294574e-06


################################################# señal de la tesis
setwd("/home/nostromo/data/pepe/")

t<-seq(0,20,0.1)

pdf("señales.pdf")
par(mfrow = c(3, 1),mar=c(4, 4, 0.8, 0.8))
y<-3*(sin(t))+sin(4*(t))
plot(t,y,type="l",, col="red", xlab="3*(seno(x))+seno(4(x))")
y<-3*(sin(t))
plot(t,y,type="l", col="blue", xlab="3*(seno(x))")
y<-sin(4*(t))
plot(t,y,type="l", col="green", xlab="seno(4(x))")
dev.off()

################################################################3
###################################################################
ANOTACION DE histona ancha



