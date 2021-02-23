source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
source("/home/nostromo/data/01_Rscripts/A_Funciones/funcionesVik.R")

library(dplyr)
library(IRanges)
library(ChIPpeakAnno)
library(GenomicRanges)


roadmap_ref<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3/hg19/E001-H3K36me3_sorted.broadPeak")

roadmap_ref<-roadmap_ref[,c(1:4)]

roadmap_ref<-setNames(roadmap_ref,c("chr","start","end","names"))

roadmap_ref$names<-paste(roadmap_ref$names)

roadmap_ranges<-GRanges(roadmap_ref)


### MACS

MACS<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3_hg19/MACS2/MACS2_h3k36me3_hg19_peaks.broadPeak")

MACS<-MACS[,c(1:4)]

MACS<-setNames(MACS,c("chr","start","end","names"))

MACS$names<-paste(MACS$names)

MACS_ranges<-GRanges(MACS)


### BCP

BCP<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3_hg19/BCP/BCP_results_h3k36me3_hg19.bed_region.bed")

BCP<-BCP[,c(1:4)]

BCP<-setNames(BCP,c("chr","start","end","names"))

BCP$names<-paste(BCP$names)

BCP_ranges<-GRanges(BCP)

### CCAT

CCAT<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3_hg19/CCAT/CCAT_results_h3k36me3_hg19.bed_region.bed")

CCAT<-CCAT[,c(1:4)]

CCAT<-setNames(CCAT,c("chr","start","end","names"))

CCAT$names<-paste(CCAT$names)

CCAT_ranges<-GRanges(CCAT)

### SICER

SICER<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3_hg19/SICER/CHIP_h3k36me3_sorted-W200-G600-islands-summary-FDR001.bed")

SICER<-SICER[,c(1:4)]

SICER<-setNames(SICER,c("chr","start","end","names"))

SICER$names<-paste(SICER$names)

SICER_ranges<-GRanges(SICER)

#### PEAKZCL

PEAKZCL<-read.table("/home/nostromo/data/pepe/39_h3k36me3_hg19/chip/OUTPUT24/PeakZCL_peaks_zCros_6_clus_10_norm_SES_min_100_max_100000.sorted_final.bed")
PEAKZCL<-read.table("/home/nostromo/data/pepe/39_h3k36me3_hg19/chip/OUTPUT25/PeakZCL_peaks_zCros_4_clus_10_norm_SES_min_100_max_100000.sorted_final.bed")
PEAKZCL<-read.table("/home/nostromo/data/pepe/39_h3k36me3_hg19/chip/OUTPUT26/PeakZCL_peaks_zCros_6_clus_10_norm_SES_min_100_max_500000.sorted_final.bed")
PEAKZCL<-read.table("/home/nostromo/data/pepe/39_h3k36me3_hg19/chip/OUTPUT26/PeakZCL_peaks_zCros_4_clus_10_norm_SES_min_100_max_10000000.sorted_final.bed")

PEAKZCL<-PEAKZCL[,c(1:4)]

PEAKZCL<-setNames(PEAKZCL,c("chr","start","end","names"))

PEAKZCL$names<-paste(PEAKZCL$names)

PEAKZCL_ranges<-GRanges(PEAKZCL)

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

roadmap_MACS<-ANNOTATION(MACS_ranges,roadmap_ranges)
#roadmap_BAYES<-ANNOTATION(BAYES_ranges,roadmap_ranges)
roadmap_SICER<-ANNOTATION(SICER_ranges,roadmap_ranges)
roadmap_BCP<-ANNOTATION(BCP_ranges,roadmap_ranges)
roadmap_CCAT<-ANNOTATION(CCAT_ranges,roadmap_ranges)
roadmap_PEAKZCL<-ANNOTATION(PEAKZCL_ranges,roadmap_ranges)


METRICS<-function(anotacion_exp,all_method,dataframe_RES){
    all_names_exp<-paste(anotacion_exp$names)
    all_names_exp<-as.numeric(lapply(strsplit(paste(all_names_exp),"Rank_"),"[",2))
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



PEAKZCL_roadmap_ref<-METRICS(roadmap_ref,PEAKZCL,roadmap_PEAKZCL)
PEAKZCL_roadmap_ref$PEAK_CALLER<-"PEAKZCL"
MACS_roadmap_ref<-METRICS(roadmap_ref,MACS,roadmap_MACS)
MACS_roadmap_ref$PEAK_CALLER<-"MACS"
BCP_roadmap_ref<-METRICS(roadmap_ref,BCP,roadmap_BCP)
BCP_roadmap_ref$PEAK_CALLER<-"BCP"
SICER_roadmap_ref<-METRICS(roadmap_ref,SICER,roadmap_SICER)
SICER_roadmap_ref$PEAK_CALLER<-"SICER"
CCAT_roadmap_ref<-METRICS(roadmap_ref,CCAT,roadmap_CCAT)
CCAT_roadmap_ref$PEAK_CALLER<-"CCAT"
results_roadmap_ref<-rbind(MACS_roadmap_ref,BCP_roadmap_ref,CCAT_roadmap_ref,SICER_roadmap_ref,PEAKZCL_roadmap_ref)
results_roadmap_ref$EXPERT<-"roadmap_ref"
write.table(results_roadmap_ref,"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/results_roadmap_ref_vs_METHOD.txt",sep="\t", col.names=TRUE, row.names=FALSE)
write.table(results_roadmap_ref,"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/results_roadmap_ref_vs_METHOD_ZCL4_min100_max_100000.txt",sep="\t", col.names=TRUE, row.names=FALSE)
write.table(results_roadmap_ref,"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/results_roadmap_ref_vs_METHOD_ZCL4_min100_max_10000000.txt",sep="\t", col.names=TRUE, row.names=FALSE)
