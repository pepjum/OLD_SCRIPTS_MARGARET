source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVik.R")

library(dplyr)
library(IRanges)
library(ChIPpeakAnno)
library(GenomicRanges)


exp1<-read.table("/home/nostromo/data/pepe/BEDS_anotaciones_manuales_chr16_h3k4me3/regiones_PEPE_h3k4me3_buena.bed")
exp2<-read.table("/home/nostromo/data/pepe/BEDS_anotaciones_manuales_chr16_h3k4me3/regions_eli_todo.bed")
exp3<-read.table("/home/nostromo/data/pepe/BEDS_anotaciones_manuales_chr16_h3k4me3/CHIP_h3k4me3_T_regions_macarena.bed")

### exp1 preparacion

exp1<-exp1[,c(1:3)]
exp1<-unique(exp1)
exp1$V4<-paste0("exp1_peak_",seq(1:nrow(exp1)))
names(exp1)<-c("chr","start","end","names")


EXP1_ranges<-GRanges(exp1)

#exp2 preparacion

exp2<-exp2[,c(1:3)]
exp2<-unique(exp2)
exp2$V4<-paste0("exp2_peak_",seq(1:nrow(exp2)))
names(exp2)<-c("chr","start","end","names")


exp2_ranges<-GRanges(exp2)

#exp3 preparacion

exp3<-exp3[,c(1:3)]
exp3<-unique(exp3)
exp3$V4<-paste0("exp3_peak_",seq(1:nrow(exp3)))
names(exp3)<-c("chr","start","end","names")

exp3_ranges<-GRanges(exp3)

#generar universo
df_universo<-rbind(exp1[,c(1:3)],exp2[,c(1:3)], exp3[,c(1:3)])
df_universo<-unique(df_universo)

universo_ranges<-GRanges(df_universo)
universo_clustered<-reduce(universo_ranges, min.gapwidth=1)
universo_clustered<-as.data.frame(universo_clustered)
universo_clustered<-universo_clustered[,c(1:3)]
universo_clustered$V4<-paste0("universo_peak_", seq(1:nrow(universo_clustered)))
names(universo_clustered)<-c("chr","start","end","names")

universo_ranges<-GRanges(universo_clustered)


ANNOTATION<-function(METHOD_RANGE,EXP_RANGE){
    peaks_overlapped_EXP_METHOD<-annotatePeakInBatch(METHOD_RANGE, EXP_RANGE, output = "both")
    peaks_from_EXP_matched_METHOD<-paste(peaks_overlapped_EXP_METHOD$feature)
    name_peak_METHOD_HIT<-paste(peaks_overlapped_EXP_METHOD$names)
    location_EXP_matched_METHOD<-peaks_overlapped_EXP_METHOD$insideFeature
    METHOD_df<-data.frame("EXP"=peaks_from_EXP_matched_METHOD,"peaks_METHOD"=name_peak_METHOD_HIT,"location"=location_EXP_matched_METHOD)
    METHOD_df_overlap<-METHOD_df[which(METHOD_df$location=="overlap"),]
    METHOD_df_overlapend<-METHOD_df[which(METHOD_df$location=="overlapEnd"),]
    METHOD_df_overlapstart<-METHOD_df[which(METHOD_df$location=="overlapStart"),]
    METHOD_df_includeFeature<-METHOD_df[which(METHOD_df$location=="includeFeature"),]
    METHOD_df_inside<-METHOD_df[which(METHOD_df$location=="inside"),]

    METHOD_df_res<-rbind(METHOD_df_overlap,METHOD_df_overlapend,METHOD_df_overlapstart,METHOD_df_includeFeature,METHOD_df_inside)
    return(METHOD_df_res)
}

EXP1_universo<-ANNOTATION(EXP1_ranges, universo_ranges)
exp2_universo<-ANNOTATION(exp2_ranges, universo_ranges)
exp3_universo<-ANNOTATION(exp3_ranges, universo_ranges)

exp1_lista<-paste(EXP1_universo$UNI)
exp2_lista<-paste(exp2_universo$UNI)
exp3_lista<-paste(exp3_universo$UNI)

pdf("/home/nostromo/data/pepe/01_RESULTADOS_CHIPSEQ_h3k4me3_ABR20/Compare3lists_EXP_universe_peak_level.pdf")
compare3List(exp1_lista,exp2_lista,exp3_lista,"EXP1","EXP2","EXP3","PEAK LEVEL")
dev.off()


library(ChIPpeakAnno)


parseENCODE <- function(x) {

	tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
	tmp <- tmp[tmp != ""]
	tmp <- tmp[seq(2, 12, 2)]
	return(tmp)
}

### AHORA CARGAMOS GENCODE y SELECCIONAMOS LAS COLUMNAS QUE NOS INTERESAN PARA ESTE TIPO DE ANALISIS
genecode<-"/home/nostromo/data/00_References/gencode_gtf/gencode.v34.annotation.gtf"

genecodev25_all <- read.table(genecode, skip = 5, header = FALSE, sep = "\t")

genecodev25 <- genecodev25_all[genecodev25_all$V3 == "gene", ]
genecodev25_tmp <- apply(as.data.frame(genecodev25[,9]), 1, parseENCODE)
genecodev25_tmp <- t(genecodev25_tmp)

colnames(genecodev25_tmp) <- c("gene_id", "gene_type", "gene_name", "gene_status", "level", "havana_gene")
genecodev25_Annot <- data.frame(genecodev25[,1:8], genecodev25_tmp[, c(1, 2,3, 4, 5,6)])

genecodev25_Annot_Gene <- genecodev25_Annot

G25AnnotChIP <- unique(genecodev25_Annot[genecodev25_Annot[, "V3"] == "gene", c("V1", "V4", "V5", "V7", "gene_id")])
colnames(G25AnnotChIP) <- c("chr", "start", "end", "strand", "gene_id")
G25AnnotChIP_RL <- RangedData(IRanges(start = G25AnnotChIP$start, end = G25AnnotChIP$end, names = paste(G25AnnotChIP$gene_id)), space = paste(G25AnnotChIP$chr), strand = paste(G25AnnotChIP$strand))

#YA TENEMOS EL OBJETO genecodev25_Annot_Gene que contiene la información que nos interesa de Gencode25 y el objeto G25AnnotChIP_RL que es una estructura de datos que contiene a qué gen corresponde cada rango de secuencias. Necesitamos los dos

#COPIA Y PEGA LA FUNCION annotChIP

annotChIP<- function(file_bed, file_annot, file_annot_gene) {

	beddata_RL <- BED2RangedData(file_bed, header=FALSE)
	beddata_Annot <- annotatePeakInBatch(beddata_RL, AnnotationData =file_annot, output = "both", multiple = F, maxgap = 0)
	beddata_Annot.df <- as.data.frame(beddata_Annot)
	beddata_Annot_G15.df <- merge(beddata_Annot.df, file_annot_gene, by.x = 7, by.y = 9)
	beddata.Filter10 <- beddata_Annot_G15.df[((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "upstream") & (abs(beddata_Annot_G15.df$distancetoFeature)<10000)) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "inside") & (abs(beddata_Annot_G15.df$distancetoFeature)<10000)) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "overlapStart")) | ((beddata_Annot_G15.df$fromOverlappingOrNearest == "NearestLocation") & (beddata_Annot_G15.df$insideFeature == "includeFeature")),]

	return(beddata.Filter10)
}

#ZCL_bed<-data_ZCL[,c("V1","V2","V3")]

#llamamos a la función annotChIP

exp1_data_annotated<-annotChIP(exp1,G25AnnotChIP_RL,genecodev25_Annot_Gene)
exp2_data_annotated<-annotChIP(exp2,G25AnnotChIP_RL,genecodev25_Annot_Gene)
exp3_data_annotated<-annotChIP(exp3,G25AnnotChIP_RL,genecodev25_Annot_Gene)


#### Escribimos la salida a un fichero de texto
#write.table(bed_data_annotated, output_NAME, quote = FALSE, col.names=T, row.names = FALSE, sep="\t")

exp1_genes<-paste(exp1_data_annotated$gene_name)
exp1_genes<-unique(exp1_genes)
exp2_genes<-paste(exp2_data_annotated$gene_name)
exp2_genes<-unique(exp2_genes)
exp3_genes<-paste(exp3_data_annotated$gene_name)
exp3_genes<-unique(exp3_genes)

pdf("/home/nostromo/data/pepe/01_RESULTADOS_CHIPSEQ_h3k4me3_ABR20/Compare3lists_genes_level.pdf")
compare3List(exp1_genes,exp2_genes,exp3_genes,"EXP1","EXP2","EXP3","GENE LEVEL")
dev.off()
