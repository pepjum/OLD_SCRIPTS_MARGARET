source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVik.R")

library(dplyr)
library(IRanges)
library(ChIPpeakAnno)
library(GenomicRanges)



ENCODE_ref<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/homo_sapiens.GRCh38.K562.H3K36me3.ccat_histone.peaks.20190329.bed")

#no tiene cromosoma Y xq proviene de una mujer

ENCODE_ref<-ENCODE_ref[,c(1:4)]
ENCODE_ref$V4<-paste0("peak_",seq(1:nrow(ENCODE_ref)))

ENCODE_ref<-setNames(ENCODE_ref,c("chr","start","end","names"))

ENCODE_ref$names<-paste(ENCODE_ref$names)
selection<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrM")

ENCODE_sel<-ENCODE_ref[which(ENCODE_ref$chr %in% selection),]


ENCODE_ranges<-GRanges(ENCODE_sel)



MACS<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/MACS2_RESULTS/MACS2_h3k36me3_peaks_maxgap1000.broadPeak")
MACS<-read.table("home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/MACS2_RESULTS_DEF/MACS2_h3k36me3_peaks.broadPeak")
MACS<-MACS[,c(1:4)]

MACS<-setNames(MACS,c("chr","start","end","names"))

MACS$names<-paste(MACS$names)

MACS_sel<-MACS[which(MACS$chr %in% selection),]


MACS_ranges<-GRanges(MACS_sel)


### BCP

BCP<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/BCP_RESULTS/BCP_h3k36me3.bed_region_winsize10000.bed")

BCP<-BCP[,c(1:4)]

BCP<-setNames(BCP,c("chr","start","end","names"))

BCP$names<-paste(BCP$names)

BCP_sel<-BCP[which(BCP$chr %in% selection),]


BCP_ranges<-GRanges(BCP_sel)

#### BayesPeak

BayesPeak<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/BayesPeak_RESULTS/Bayespeak_h3k36me3.txt")

BayesPeak<-BayesPeak[,c(1:4)]

BayesPeak<-setNames(BayesPeak,c("chr","start","end","names"))

BayesPeak$names<-paste(BayesPeak$names)

BayesPeak_sel<-BayesPeak[which(BayesPeak$chr %in% selection),]

BayesPeak_ranges<-GRanges(BayesPeak_sel)


#CCAT

CCAT<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/CCAT_RESULTS/CCAT_h3k36me3.bed_region_winsize10000.bed")

CCAT<-CCAT[,c(1:4)]

CCAT<-setNames(CCAT,c("chr","start","end","names"))

CCAT$names<-paste(CCAT$names)

CCAT_sel<-CCAT[which(CCAT$chr %in% selection),]


CCAT_ranges<-GRanges(CCAT_sel)

### SICER

SICER<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/SICER_RESULTS/CHIP_sorted_h3k36me3-W200-G600-islands-summary-FDR0.01.bed")

SICER<-SICER[,c(1:4)]

SICER<-setNames(SICER,c("chr","start","end","names"))

SICER$names<-paste(SICER$names)

SICER_sel<-SICER[which(SICER$chr %in% selection),]


SICER_ranges<-GRanges(SICER_sel)

#### JAMM

JAMM<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/JAMM_RESULTS/peaks/selected_JAMM_h3k36me3_hg38.narrowPeak")

JAMM<-JAMM[,c(1:4)]

JAMM<-setNames(JAMM,c("chr","start","end","names"))

JAMM$names<-paste(JAMM$names)

JAMM_sel<-JAMM[which(JAMM$chr %in% selection),]

JAMM_ranges<-GRanges(JAMM_sel)

### PEAKZCL

PEAKZCL<-read.table("/home/nostromo/data/pepe/40_h3k36me3_hg38/chip/OUTPUT001/PeakZCL_peaks_zCros_4_clus_10_norm_SES_min_500_max_10000000.sorted_final.bed")

PEAKZCL<-PEAKZCL[,c(1:4)]

PEAKZCL<-setNames(PEAKZCL,c("chr","start","end","names"))

PEAKZCL$names<-paste(PEAKZCL$names)

PEAKZCL_sel<-PEAKZCL[which(PEAKZCL$chr %in% selection),]

PEAKZCL_ranges<-GRanges(PEAKZCL_sel)


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

ENCODE_MACS<-ANNOTATION(MACS_ranges,ENCODE_ranges)
ENCODE_JAMM<-ANNOTATION(JAMM_ranges,ENCODE_ranges)
ENCODE_SICER<-ANNOTATION(SICER_ranges,ENCODE_ranges)
ENCODE_BCP<-ANNOTATION(BCP_ranges,ENCODE_ranges)
ENCODE_CCAT<-ANNOTATION(CCAT_ranges,ENCODE_ranges)
ENCODE_PEAKZCL<-ANNOTATION(PEAKZCL_ranges,ENCODE_ranges)
ENCODE_BayesPeak<-ANNOTATION(BayesPeak_ranges,ENCODE_ranges)

METRICS<-function(anotacion_exp,all_method,dataframe_RES){
    all_names_exp<-paste(anotacion_exp$names)
    all_names_exp<-as.numeric(lapply(strsplit(paste(all_names_exp),"peak_"),"[",2))
    peaks_exp_MATCHED_method<-unique(as.numeric(paste(dataframe_RES[,1])))
    TRUE_POSITIVES<-length(intersect(all_names_exp,peaks_exp_MATCHED_method))
    FALSE_POSITIVES<-length(setdiff(paste(all_method[,4]),unique(paste(dataframe_RES[,2]))))
    FALSE_NEGATIVES<-length(setdiff(all_names_exp,paste(peaks_exp_MATCHED_method)))
    TPR<-TRUE_POSITIVES/(TRUE_POSITIVES + FALSE_NEGATIVES)
    PPV<-TRUE_POSITIVES/(TRUE_POSITIVES + FALSE_POSITIVES)
    FNR<-1-TPR
    FDR<-1-PPV
    F1_score<-2*(PPV*TPR/(PPV+TPR))
    results<-data.frame("TRUE_POSITIVES"=as.numeric(TRUE_POSITIVES), "FALSE_POSITIVES"=as.numeric(FALSE_POSITIVES),"FALSE_NEGATIVES"=as.numeric(FALSE_NEGATIVES),"TPR"=as.numeric(TPR),"PPV"=as.numeric(PPV),"FNR"=as.numeric(FNR),"FDR"=as.numeric(FDR),"F1_score"=as.numeric(F1_score))
    return(results)
}


PEAKZCL_ENCODE_sel<-METRICS(ENCODE_sel,PEAKZCL_sel,ENCODE_PEAKZCL)
PEAKZCL_ENCODE_sel$PEAK_CALLER<-"PEAKZCL"
MACS_ENCODE_sel<-METRICS(ENCODE_sel,MACS_sel,ENCODE_MACS)
MACS_ENCODE_sel$PEAK_CALLER<-"MACS"
BCP_ENCODE_sel<-METRICS(ENCODE_sel,BCP_sel,ENCODE_BCP)
BCP_ENCODE_sel$PEAK_CALLER<-"BCP"
SICER_ENCODE_sel<-METRICS(ENCODE_sel,SICER_sel,ENCODE_SICER)
SICER_ENCODE_sel$PEAK_CALLER<-"SICER"
CCAT_ENCODE_sel<-METRICS(ENCODE_sel,CCAT_sel,ENCODE_CCAT)
CCAT_ENCODE_sel$PEAK_CALLER<-"CCAT"
JAMM_ENCODE_sel<-METRICS(ENCODE_sel,JAMM_sel,ENCODE_JAMM)
JAMM_ENCODE_sel$PEAK_CALLER<-"JAMM"
BayesPeak_ENCODE_sel<-METRICS(ENCODE_sel,BayesPeak_sel,ENCODE_BayesPeak)
BayesPeak_ENCODE_sel$PEAK_CALLER<-"BayesPeak"

results_ENCODE_sel<-rbind(MACS_ENCODE_sel,BCP_ENCODE_sel,CCAT_ENCODE_sel,SICER_ENCODE_sel,JAMM_ENCODE_sel,BayesPeak_ENCODE_sel,PEAKZCL_ENCODE_sel)
results_ENCODE_sel$EXPERT<-"ENCODE"
write.table(results_ENCODE_sel,file="/home/nostromo/data/pepe/01_RESULTADOS_CHIPSEQ_h3k4me3_ABR20/results_h3k36me3_hg38.txt", col.names=T,quote=F , sep="\t")
write.table(results_ENCODE_sel,file="/home/nostromo/data/pepe/01_RESULTADOS_CHIPSEQ_h3k4me3_ABR20/results_h3k36me3_hg38_MACS_defecto.txt", col.names=T,quote=F , sep="\t")
