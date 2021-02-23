#### CALCULAR LAS AREAS DE LOS METODS Y LAS AMPLITUDES. CREAR GRAFICOS  IGUAL al script GENERATE_PLOTS_PAPER_ZCL_JUN19.R


### files h3k4me3


MACS_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/MACS2_RESULTS/MACS2_h3k4me3_peaks.narrowPeak"
JAMM_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/JAMM_RESULTS/peaks/selected_JAMM_h3k4me3.narrowPeak"
SICER_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/SICER_RESULTS/CHIP_sorted_h3k4me3-W200-G200-islands-summary-FDR0.01.bed"
BAYES_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/BayesPeak_RESULTS/BayesPeak_output_h3k4me3.bed"
RANGER_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/RANGER_RESULTS/RANGER_h3k4me3.bed_region.bed"
ZCL_file<-"/home/nostromo/data/pepe/36_h3k4me3/chip/OUTPUT23/PeakZCL_peaks_zCros_22_clus_10_norm_SES_min_150_max_6500.sorted_final.bed"


bedgraph_chr16_file<-"/home/nostromo/data/pepe/36_h3k4me3/chip/TMP/chip/chr16.chrom"

MACS<-read.table(MACS_file, header=F)
MACS_chr16<-MACS[which(MACS$V1=="chr16"),]
MACS_chr16$amplitud<-as.numeric(paste(MACS_chr16$V3))-as.numeric(paste(MACS_chr16$V2))

JAMM<-read.table(JAMM_file, header=F)
JAMM_chr16<-JAMM[which(JAMM$V1=="chr16"),]
JAMM_chr16$amplitud<-as.numeric(paste(JAMM_chr16$V3))-as.numeric(paste(JAMM_chr16$V2))

BAYES<-read.table(BAYES_file, header=F)
BAYES_chr16<-BAYES[which(BAYES$V1=="chr16"),]
BAYES_chr16$amplitud<-as.numeric(paste(BAYES_chr16$V3))-as.numeric(paste(BAYES_chr16$V2))

RANGER<-read.table(RANGER_file, header=F)
RANGER_chr16<-RANGER[which(RANGER$V1=="chr16"),]
RANGER_chr16$amplitud<-as.numeric(paste(RANGER_chr16$V3))-as.numeric(paste(RANGER_chr16$V2))

SICER<-read.table(SICER_file, header=F)
SICER_chr16<-SICER[which(SICER$V1=="chr16"),]
SICER_chr16$amplitud<-as.numeric(paste(SICER_chr16$V3))-as.numeric(paste(SICER_chr16$V2))

ZCL<-read.table(ZCL_file, header=F)
ZCL_chr16<-ZCL[which(ZCL$V1=="chr16"),]
ZCL_chr16$amplitud<-as.numeric(paste(ZCL_chr16$V3))-as.numeric(paste(ZCL_chr16$V2))

exp1<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/regiones_PEPE_h3k4me3_buena.bed")
exp2<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/regions_eli_todo.bed")
exp3<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/CHIP_h3k4me3_T_regions_macarena.bed")

exp1<-unique(exp1)
names(exp1)<-c("V1","V2","V3","V4")
exp2<-unique(exp2)
exp3<-unique(exp3)
exp1$amplitud<-as.numeric(exp1$V3)-as.numeric(exp1$V2)
exp2$amplitud<-as.numeric(exp2$V3)-as.numeric(exp2$V2)
exp3$amplitud<-as.numeric(exp3$V3)-as.numeric(exp3$V2)



require(ggplot2)
setwd("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/")

####    AREAS
require(data.table)

bedgraph_chr16<-fread(bedgraph_chr16_file)
bedgraph_chr16<-as.data.frame(bedgraph_chr16)


quantification<-function(dataframe,signal){
  quantified_peaks<-list()
  for (i in 1:nrow(dataframe)){
    peaks_quantified<-sum(signal$V3[dataframe$V2[i]:dataframe$V3[i]])
    quantified_peaks[[i]]<-peaks_quantified
  }
  return(quantified_peaks)
}

MACS2_areas<-quantification(MACS_chr16,bedgraph_chr16)
MACS_chr16$areas<-paste(MACS2_areas)
MACS_chr16$log_areas<-log2(as.numeric(paste(MACS_chr16$areas)))

JAMM_areas<-quantification(JAMM_chr16,bedgraph_chr16)
JAMM_chr16$areas<-paste(JAMM_areas)
JAMM_chr16$log_areas<-log2(as.numeric(paste(JAMM_chr16$areas)))

BAYES_areas<-quantification(BAYES_chr16,bedgraph_chr16)
BAYES_chr16$areas<-paste(BAYES_areas)
BAYES_chr16$log_areas<-log2(as.numeric(paste(BAYES_chr16$areas)))

RANGER_areas<-quantification(RANGER_chr16,bedgraph_chr16)
RANGER_chr16$areas<-paste(RANGER_areas)
RANGER_chr16$log_areas<-log2(as.numeric(paste(RANGER_chr16$areas)))


SICER_areas<-quantification(SICER_chr16,bedgraph_chr16)
SICER_chr16$areas<-paste(SICER_areas)
SICER_chr16$log_areas<-log2(as.numeric(paste(SICER_chr16$areas)))


ZCL_areas<-quantification(ZCL_chr16,bedgraph_chr16)
ZCL_chr16$areas<-paste(ZCL_areas)
ZCL_chr16$log_areas<-log2(as.numeric(paste(ZCL_chr16$areas)))


exp1_areas<-quantification(exp1, bedgraph_chr16)
exp1$areas<-paste(exp1_areas)
exp1$log2_areas<-log2(as.numeric(paste(exp1$areas)))

exp2_areas<-quantification(exp2, bedgraph_chr16)
exp2$areas<-paste(exp2_areas)
exp2$log2_areas<-log2(as.numeric(paste(exp2$areas)))

exp3_areas<-quantification(exp3, bedgraph_chr16)
exp3$areas<-paste(exp3_areas)
exp3$log2_areas<-log2(as.numeric(paste(exp3$areas)))

write.table(exp1, file="exp1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(exp2, file="exp2_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(exp3, file="exp3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(ZCL_chr16, file="PeakZCL_chr16__h3k4me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(MACS_chr16, file="MACS_chr16_h3k4me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(SICER_chr16, file="SICER_chr16_h3k4me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(RANGER_chr16, file="RANGER_chr16_h3k4me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(JAMM_chr16, file="JAMM_chr16_h3k4me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(BAYES_chr16, file="BayesPeak_chr16_h3k4me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")



###########################################################
#########################################################

#EGR1

MACS_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/MACS2_RESULTS/MACS2_EGR1_peaks.narrowPeak"
JAMM_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/JAMM_RESULTS/peaks/selected_JAMM_EGR1.narrowPeak"
SICER_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/SICER_RESULTS/CHIP_sorted_EGR1-W200-G200-islands-summary-FDR0.01.bed"
BAYES_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/BayesPeak_RESULTS/BayesPeak_output_EGR1.bed"
RANGER_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/RANGER_RESULTS/RANGER_EGR1.bed_region.bed"
ZCL_file<-"/home/nostromo/data/pepe/35_CHIPSEQ_EGR1_BAMS/chip/OUTPUT24/PeakZCL_peaks_zCros_22_clus_1_norm_SES_min_100_max_1500.sorted_final.bed"


setwd("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/")


MACS<-read.table(MACS_file, header=F)
MACS_chr16<-MACS[which(MACS$V1=="chr16"),]
MACS_chr16$amplitud<-as.numeric(paste(MACS_chr16$V3))-as.numeric(paste(MACS_chr16$V2))

JAMM<-read.table(JAMM_file, header=F)
JAMM_chr16<-JAMM[which(JAMM$V1=="chr16"),]
JAMM_chr16$amplitud<-as.numeric(paste(JAMM_chr16$V3))-as.numeric(paste(JAMM_chr16$V2))

BAYES<-read.table(BAYES_file, header=F)
BAYES_chr16<-BAYES[which(BAYES$V1=="chr16"),]
BAYES_chr16$amplitud<-as.numeric(paste(BAYES_chr16$V3))-as.numeric(paste(BAYES_chr16$V2))

RANGER<-read.table(RANGER_file, header=F)
RANGER_chr16<-RANGER[which(RANGER$V1=="chr16"),]
RANGER_chr16$amplitud<-as.numeric(paste(RANGER_chr16$V3))-as.numeric(paste(RANGER_chr16$V2))

SICER<-read.table(SICER_file, header=F)
SICER_chr16<-SICER[which(SICER$V1=="chr16"),]
SICER_chr16$amplitud<-as.numeric(paste(SICER_chr16$V3))-as.numeric(paste(SICER_chr16$V2))

ZCL<-read.table(ZCL_file, header=F)
ZCL_chr16<-ZCL[which(ZCL$V1=="chr16"),]
ZCL_chr16$amplitud<-as.numeric(paste(ZCL_chr16$V3))-as.numeric(paste(ZCL_chr16$V2))


bedgraph_chr16_file<-"/home/nostromo/data/pepe/35_CHIPSEQ_EGR1_BAMS/chip/TMP/chip/chr16.chrom"
bedgraph_chr16<-fread(bedgraph_chr16_file)
bedgraph_chr16<-as.data.frame(bedgraph_chr16)


MACS2_areas<-quantification(MACS_chr16,bedgraph_chr16)
MACS_chr16$areas<-paste(MACS2_areas)
MACS_chr16$log_areas<-log2(as.numeric(paste(MACS_chr16$areas)))

JAMM_areas<-quantification(JAMM_chr16,bedgraph_chr16)
JAMM_chr16$areas<-paste(JAMM_areas)
JAMM_chr16$log_areas<-log2(as.numeric(paste(JAMM_chr16$areas)))

BAYES_areas<-quantification(BAYES_chr16,bedgraph_chr16)
BAYES_chr16$areas<-paste(BAYES_areas)
BAYES_chr16$log_areas<-log2(as.numeric(paste(BAYES_chr16$areas)))

RANGER_areas<-quantification(RANGER_chr16,bedgraph_chr16)
RANGER_chr16$areas<-paste(RANGER_areas)
RANGER_chr16$log_areas<-log2(as.numeric(paste(RANGER_chr16$areas)))

SICER_areas<-quantification(SICER_chr16,bedgraph_chr16)
SICER_chr16$areas<-paste(SICER_areas)
SICER_chr16$log_areas<-log2(as.numeric(paste(SICER_chr16$areas)))

ZCL_areas<-quantification(ZCL_chr16,bedgraph_chr16)
ZCL_chr16$areas<-paste(ZCL_areas)
ZCL_chr16$log_areas<-log2(as.numeric(paste(ZCL_chr16$areas)))

write.table(ZCL_chr16, file="PeakZCL_chr16__EGR1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(MACS_chr16, file="MACS_chr16_EGR1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(SICER_chr16, file="SICER_chr16_EGR1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(RANGER_chr16, file="RANGER_chr16_EGR1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(JAMM_chr16, file="JAMM_chr16_EGR1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(BAYES_chr16, file="BayesPeak_chr16_EGR1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")


###########################################################
#########################################################

#JUND

MACS_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/MACS2_RESULTS/MACS2_JUND_peaks.narrowPeak"
JAMM_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/JAMM_RESULTS/peaks/selected_JAMM_JUND.narrowPeak"
SICER_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/SICER_RESULTS/CHIP_sorted_JUND-W200-G200-islands-summary-FDR0.01.bed"
BAYES_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/BayesPeak_RESULTS/BayesPeak_output_JUND.bed"
RANGER_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/RANGER_RESULTS/RANGER_JUND.bed_region.bed"
ZCL_file<-"/home/nostromo/data/pepe/37_JUND/chip/OUTPUT24/PeakZCL_peaks_zCros_24_clus_1_norm_SES_min_60_max_1000.sorted_final.bed"


setwd("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/")


MACS<-read.table(MACS_file, header=F)
MACS_chr16<-MACS[which(MACS$V1=="chr16"),]
MACS_chr16$amplitud<-as.numeric(paste(MACS_chr16$V3))-as.numeric(paste(MACS_chr16$V2))

JAMM<-read.table(JAMM_file, header=F)
JAMM_chr16<-JAMM[which(JAMM$V1=="chr16"),]
JAMM_chr16$amplitud<-as.numeric(paste(JAMM_chr16$V3))-as.numeric(paste(JAMM_chr16$V2))

BAYES<-read.table(BAYES_file, header=F)
BAYES_chr16<-BAYES[which(BAYES$V1=="chr16"),]
BAYES_chr16$amplitud<-as.numeric(paste(BAYES_chr16$V3))-as.numeric(paste(BAYES_chr16$V2))

RANGER<-read.table(RANGER_file, header=F)
RANGER_chr16<-RANGER[which(RANGER$V1=="chr16"),]
RANGER_chr16$amplitud<-as.numeric(paste(RANGER_chr16$V3))-as.numeric(paste(RANGER_chr16$V2))

SICER<-read.table(SICER_file, header=F)
SICER_chr16<-SICER[which(SICER$V1=="chr16"),]
SICER_chr16$amplitud<-as.numeric(paste(SICER_chr16$V3))-as.numeric(paste(SICER_chr16$V2))

ZCL<-read.table(ZCL_file, header=F)
ZCL_chr16<-ZCL[which(ZCL$V1=="chr16"),]
ZCL_chr16$amplitud<-as.numeric(paste(ZCL_chr16$V3))-as.numeric(paste(ZCL_chr16$V2))


bedgraph_chr16_file<-"/home/nostromo/data/pepe/37_JUND/chip/TMP/chip/chr16.chrom"
bedgraph_chr16<-fread(bedgraph_chr16_file)
bedgraph_chr16<-as.data.frame(bedgraph_chr16)


MACS2_areas<-quantification(MACS_chr16,bedgraph_chr16)
MACS_chr16$areas<-paste(MACS2_areas)
MACS_chr16$log_areas<-log2(as.numeric(paste(MACS_chr16$areas)))

JAMM_areas<-quantification(JAMM_chr16,bedgraph_chr16)
JAMM_chr16$areas<-paste(JAMM_areas)
JAMM_chr16$log_areas<-log2(as.numeric(paste(JAMM_chr16$areas)))

BAYES_areas<-quantification(BAYES_chr16,bedgraph_chr16)
BAYES_chr16$areas<-paste(BAYES_areas)
BAYES_chr16$log_areas<-log2(as.numeric(paste(BAYES_chr16$areas)))

RANGER_areas<-quantification(RANGER_chr16,bedgraph_chr16)
RANGER_chr16$areas<-paste(RANGER_areas)
RANGER_chr16$log_areas<-log2(as.numeric(paste(RANGER_chr16$areas)))

SICER_areas<-quantification(SICER_chr16,bedgraph_chr16)
SICER_chr16$areas<-paste(SICER_areas)
SICER_chr16$log_areas<-log2(as.numeric(paste(SICER_chr16$areas)))

ZCL_areas<-quantification(ZCL_chr16,bedgraph_chr16)
ZCL_chr16$areas<-paste(ZCL_areas)
ZCL_chr16$log_areas<-log2(as.numeric(paste(ZCL_chr16$areas)))

write.table(ZCL_chr16, file="PeakZCL_chr16__JUND_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(MACS_chr16, file="MACS_chr16_JUND_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(SICER_chr16, file="SICER_chr16_JUND_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(RANGER_chr16, file="RANGER_chr16_JUND_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(JAMM_chr16, file="JAMM_chr16_JUND_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(BAYES_chr16, file="BayesPeak_chr16_JUND_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")



###########################################################
#########################################################

#SP1

MACS_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/MACS2_RESULTS/MACS2_SP1_peaks.narrowPeak"
JAMM_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/JAMM_RESULTS/peaks/selected_JAMM_SP1.narrowPeak"
SICER_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/SICER_RESULTS/CHIP_sorted_SP1-W200-G200-islands-summary-FDR0.01.bed"
BAYES_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/BAYESPEAK_RESULTS/BayesPeak_output_SP1.bed"
RANGER_file<-"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/RANGER_RESULTS/RANGER_SP1.bed_region.bed"
ZCL_file<-"/home/nostromo/data/pepe/38_SP1/chip/OUTPUT23/PeakZCL_peaks_zCros_24_clus_1_norm_SES_min_120_max_1000.sorted_final.bed"


setwd("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/")


MACS<-read.table(MACS_file, header=F)
MACS_chr16<-MACS[which(MACS$V1=="chr16"),]
MACS_chr16$amplitud<-as.numeric(paste(MACS_chr16$V3))-as.numeric(paste(MACS_chr16$V2))

JAMM<-read.table(JAMM_file, header=F)
JAMM_chr16<-JAMM[which(JAMM$V1=="chr16"),]
JAMM_chr16$amplitud<-as.numeric(paste(JAMM_chr16$V3))-as.numeric(paste(JAMM_chr16$V2))

BAYES<-read.table(BAYES_file, header=F)
BAYES_chr16<-BAYES[which(BAYES$V1=="chr16"),]
BAYES_chr16$amplitud<-as.numeric(paste(BAYES_chr16$V3))-as.numeric(paste(BAYES_chr16$V2))

RANGER<-read.table(RANGER_file, header=F)
RANGER_chr16<-RANGER[which(RANGER$V1=="chr16"),]
RANGER_chr16$amplitud<-as.numeric(paste(RANGER_chr16$V3))-as.numeric(paste(RANGER_chr16$V2))

SICER<-read.table(SICER_file, header=F)
SICER_chr16<-SICER[which(SICER$V1=="chr16"),]
SICER_chr16$amplitud<-as.numeric(paste(SICER_chr16$V3))-as.numeric(paste(SICER_chr16$V2))

ZCL<-read.table(ZCL_file, header=F)
ZCL_chr16<-ZCL[which(ZCL$V1=="chr16"),]
ZCL_chr16$amplitud<-as.numeric(paste(ZCL_chr16$V3))-as.numeric(paste(ZCL_chr16$V2))


bedgraph_chr16_file<-"/home/nostromo/data/pepe/38_SP1/chip/TMP/chip/chr16.chrom"
bedgraph_chr16<-fread(bedgraph_chr16_file)
bedgraph_chr16<-as.data.frame(bedgraph_chr16)


MACS2_areas<-quantification(MACS_chr16,bedgraph_chr16)
MACS_chr16$areas<-paste(MACS2_areas)
MACS_chr16$log_areas<-log2(as.numeric(paste(MACS_chr16$areas)))

JAMM_areas<-quantification(JAMM_chr16,bedgraph_chr16)
JAMM_chr16$areas<-paste(JAMM_areas)
JAMM_chr16$log_areas<-log2(as.numeric(paste(JAMM_chr16$areas)))

BAYES_areas<-quantification(BAYES_chr16,bedgraph_chr16)
BAYES_chr16$areas<-paste(BAYES_areas)
BAYES_chr16$log_areas<-log2(as.numeric(paste(BAYES_chr16$areas)))

RANGER_areas<-quantification(RANGER_chr16,bedgraph_chr16)
RANGER_chr16$areas<-paste(RANGER_areas)
RANGER_chr16$log_areas<-log2(as.numeric(paste(RANGER_chr16$areas)))

SICER_areas<-quantification(SICER_chr16,bedgraph_chr16)
SICER_chr16$areas<-paste(SICER_areas)
SICER_chr16$log_areas<-log2(as.numeric(paste(SICER_chr16$areas)))

ZCL_areas<-quantification(ZCL_chr16,bedgraph_chr16)
ZCL_chr16$areas<-paste(ZCL_areas)
ZCL_chr16$log_areas<-log2(as.numeric(paste(ZCL_chr16$areas)))

write.table(ZCL_chr16, file="PeakZCL_chr16__SP1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(MACS_chr16, file="MACS_chr16_SP1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(SICER_chr16, file="SICER_chr16_SP1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(RANGER_chr16, file="RANGER_chr16_SP1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(JAMM_chr16, file="JAMM_chr16_SP1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(BAYES_chr16, file="BayesPeak_chr16_SP1_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")




###########################################################
#########################################################

#h3k36me3

MACS_file<-"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/MACS2/MACS2_h3k36me3_peaks.broadPeak"
JAMM_file<-"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/JAMM/peaks/selected_JAMM_h3k36me3_NOV20_hg19.narrowPeak"
SICER_file<-"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/SICER/CHIP_h3k36me3_hg19_NOV20_sorted-W200-G600-islands-summary-FDR0.1.bed"
BAYES_file<-"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/BayesPeak_h3k36me3_hg19_NOV20.bed"
BCP_file<-"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/BCP_0.01/BCP_0.01_h3k36me3_hg19_NOV.bed_region.bed"
CCAT_file<-"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/CCAT/CCAT_h3k36me3_hg19_NOV20.bed_region.bed"
ZCL_file<-"/home/nostromo/data/pepe/39_h3k36me3_hg19/chip/OUTPUT26/PeakZCL_peaks_zCros_4_clus_10_norm_SES_min_100_max_10000000.sorted_final.bed"


setwd("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/")


MACS<-read.table(MACS_file, header=F)
MACS_chr16<-MACS[which(MACS$V1=="chr16"),]
MACS_chr16$amplitud<-as.numeric(paste(MACS_chr16$V3))-as.numeric(paste(MACS_chr16$V2))

JAMM<-read.table(JAMM_file, header=F)
JAMM_chr16<-JAMM[which(JAMM$V1=="chr16"),]
JAMM_chr16$amplitud<-as.numeric(paste(JAMM_chr16$V3))-as.numeric(paste(JAMM_chr16$V2))

BAYES<-read.table(BAYES_file, header=F)
BAYES_chr16<-BAYES[which(BAYES$V1=="chr16"),]
BAYES_chr16$amplitud<-as.numeric(paste(BAYES_chr16$V3))-as.numeric(paste(BAYES_chr16$V2))

BCP<-read.table(BCP_file, header=F)
BCP_chr16<-BCP[which(BCP$V1=="chr16"),]
BCP_chr16$amplitud<-as.numeric(paste(BCP_chr16$V3))-as.numeric(paste(BCP_chr16$V2))

CCAT<-read.table(CCAT_file, header=F)
CCAT_chr16<-CCAT[which(CCAT$V1=="chr16"),]
CCAT_chr16$amplitud<-as.numeric(paste(CCAT_chr16$V3))-as.numeric(paste(CCAT_chr16$V2))


SICER<-read.table(SICER_file, header=F)
SICER_chr16<-SICER[which(SICER$V1=="chr16"),]
SICER_chr16$amplitud<-as.numeric(paste(SICER_chr16$V3))-as.numeric(paste(SICER_chr16$V2))

ZCL<-read.table(ZCL_file, header=F)
ZCL_chr16<-ZCL[which(ZCL$V1=="chr16"),]
ZCL_chr16$amplitud<-as.numeric(paste(ZCL_chr16$V3))-as.numeric(paste(ZCL_chr16$V2))


bedgraph_chr16_file<-"/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/chip_bam/TMP/chip/chr16.chrom"
bedgraph_chr16<-fread(bedgraph_chr16_file)
bedgraph_chr16<-as.data.frame(bedgraph_chr16)


MACS2_areas<-quantification(MACS_chr16,bedgraph_chr16)
MACS_chr16$areas<-paste(MACS2_areas)
MACS_chr16$log_areas<-log2(as.numeric(paste(MACS_chr16$areas)))

JAMM_areas<-quantification(JAMM_chr16,bedgraph_chr16)
JAMM_chr16$areas<-paste(JAMM_areas)
JAMM_chr16$log_areas<-log2(as.numeric(paste(JAMM_chr16$areas)))

BAYES_areas<-quantification(BAYES_chr16,bedgraph_chr16)
BAYES_chr16$areas<-paste(BAYES_areas)
BAYES_chr16$log_areas<-log2(as.numeric(paste(BAYES_chr16$areas)))

BCP_areas<-quantification(BCP_chr16,bedgraph_chr16)
BCP_chr16$areas<-paste(BCP_areas)
BCP_chr16$log_areas<-log2(as.numeric(paste(BCP_chr16$areas)))

CCAT_areas<-quantification(CCAT_chr16,bedgraph_chr16)
CCAT_chr16$areas<-paste(CCAT_areas)
CCAT_chr16$log_areas<-log2(as.numeric(paste(CCAT_chr16$areas)))


SICER_areas<-quantification(SICER_chr16,bedgraph_chr16)
SICER_chr16$areas<-paste(SICER_areas)
SICER_chr16$log_areas<-log2(as.numeric(paste(SICER_chr16$areas)))

ZCL_areas<-quantification(ZCL_chr16,bedgraph_chr16)
ZCL_chr16$areas<-paste(ZCL_areas)
ZCL_chr16$log_areas<-log2(as.numeric(paste(ZCL_chr16$areas)))

write.table(ZCL_chr16, file="PeakZCL_chr16__h3k36me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(MACS_chr16, file="MACS_chr16_h3k36me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(SICER_chr16, file="SICER_chr16_h3k36me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(BCP_chr16, file="BCP_chr16_h3k36me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(CCAT_chr16, file="CCAT_chr16_h3k36me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(JAMM_chr16, file="JAMM_chr16_h3k36me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")
write.table(BAYES_chr16, file="BayesPeak_chr16_h3k36me3_areas_anchuras.bed", col.names=T, row.names=F, quote=F, sep="\t")









