#### CALCULAR LAS AREAS DE LOS METODS Y LAS AMPLITUDES. CREAR GRAFICOS  IGUAL al script GENERATE_PLOTS_PAPER_ZCL_JUN19.R


### files h3k4me3


MACS_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/MACS2/CHIP_h3k4me3_peaks.narrowPeak"
JAMM_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/JAMM/selected_JAMM.narrowPeak"
BAYES_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/BAYES/BayesPeak.bed"
CCAT_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/CCAT/CCAT_h3k4me3.bed_region_passed.bed"
BCP_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/BCP/h3k4me3_bcp.bed_region.bed"
SICER_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/SICER/CHIP_h3k4me3-W200-G200-FDR0.01-island.bed"
ZCL_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/ZCL/ZCL_peaks.bed"

bedgraph_chr16_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/cromosomas_h3k4me3/chr16_h3k4me3.chrom"

MACS<-read.table(MACS_file, header=F)
MACS_chr16<-MACS[which(MACS$V1=="chr16"),]
MACS_chr16$amplitud<-as.numeric(paste(MACS_chr16$V3))-as.numeric(paste(MACS_chr16$V2))

JAMM<-read.table(JAMM_file, header=F)
JAMM_chr16<-JAMM[which(JAMM$V1=="chr16"),]
JAMM_chr16$amplitud<-as.numeric(paste(JAMM_chr16$V3))-as.numeric(paste(JAMM_chr16$V2))

BAYES<-read.table(BAYES_file, header=F)
BAYES_chr16<-BAYES[which(BAYES$V1=="chr16"),]
BAYES_chr16$amplitud<-as.numeric(paste(BAYES_chr16$V3))-as.numeric(paste(BAYES_chr16$V2))

CCAT<-read.table(CCAT_file, header=F)
CCAT_chr16<-CCAT[which(CCAT$V1=="chr16"),]
CCAT_chr16$amplitud<-as.numeric(paste(CCAT_chr16$V3))-as.numeric(paste(CCAT_chr16$V2))

BCP<-read.table(BCP_file, header=F)
BCP_chr16<-BCP[which(BCP$V1=="chr16"),]
BCP_chr16$amplitud<-as.numeric(paste(BCP_chr16$V3))-as.numeric(paste(BCP_chr16$V2))

SICER<-read.table(SICER_file, header=F)
SICER_chr16<-SICER[which(SICER$V1=="chr16"),]
SICER_chr16$amplitud<-as.numeric(paste(SICER_chr16$V3))-as.numeric(paste(SICER_chr16$V2))

ZCL<-read.table(ZCL_file, header=F)
ZCL_chr16<-ZCL[which(ZCL$V1=="chr16"),]
ZCL_chr16$amplitud<-as.numeric(paste(ZCL_chr16$V3))-as.numeric(paste(ZCL_chr16$V2))

exp1<-read.table("/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/exp1_w_names.bed")
exp2<-read.table("/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/exp2_w_names.bed")
exp3<-read.table("/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/CHIP_h3k4me3_T_regions_macarena.bed")

exp1<-unique(exp1)
exp2<-unique(exp2)
exp3<-unique(exp3)
exp1$amplitud<-as.numeric(exp1$V3)-as.numeric(exp1$V2)
exp2$amplitud<-as.numeric(exp2$V3)-as.numeric(exp2$V2)
exp3$amplitud<-as.numeric(exp3$V3)-as.numeric(exp3$V2)



require(ggplot2)

plotter <- data.frame('amplitude' = NULL, 'Method' = NULL)

plotter <- rbind(plotter, data.frame('amplitude' = MACS_chr16$amplitud, 'Method' = rep('MACS2', nrow(MACS_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = JAMM_chr16$amplitud, 'Method' = rep('JAMM', nrow(JAMM_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = BAYES_chr16$amplitud, 'Method' = rep('BayesPeak', nrow(BAYES_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = CCAT_chr16$amplitud, 'Method' = rep('CCAT', nrow(CCAT_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = BCP_chr16$amplitud, 'Method' = rep('BCP', nrow(BCP_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = SICER_chr16$amplitud, 'Method' = rep('SICER', nrow(SICER_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = ZCL_chr16$amplitud, 'Method' = rep('PeakZCL', nrow(ZCL_chr16))))

pdf("violin_plot_amplitudes_metodos.pdf")
ggplot(plotter, aes(x=Method, y=amplitude)) + geom_violin() +geom_boxplot(width=0.1, outlier.shape = NA) +  theme_bw()
dev.off()

plotter <- rbind(plotter, data.frame('amplitude' = exp1$amplitud, 'Method' = rep('EXP1', nrow(exp1))))
plotter <- rbind(plotter, data.frame('amplitude' = exp2$amplitud, 'Method' = rep('EXP2', nrow(exp2))))
plotter <- rbind(plotter, data.frame('amplitude' = exp3$amplitud, 'Method' = rep('EXP3', nrow(exp3))))

pdf("violin_plot_amplitudes_ALL.pdf")
ggplot(plotter, aes(x=Method, y=amplitude)) + geom_violin() +geom_boxplot(width=0.1, outlier.shape = NA) +  theme_bw()
dev.off()


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
MACS_chr16$log_areas<-log10(as.numeric(paste(MACS_chr16$areas)))

JAMM_areas<-quantification(JAMM_chr16,bedgraph_chr16)
JAMM_chr16$areas<-paste(JAMM_areas)
JAMM_chr16$log_areas<-log10(as.numeric(paste(JAMM_chr16$areas)))

BAYES_areas<-quantification(BAYES_chr16,bedgraph_chr16)
BAYES_chr16$areas<-paste(BAYES_areas)
BAYES_chr16$log_areas<-log10(as.numeric(paste(BAYES_chr16$areas)))


CCAT_areas<-quantification(CCAT_chr16,bedgraph_chr16)
CCAT_chr16$areas<-paste(CCAT_areas)
CCAT_chr16$log_areas<-log10(as.numeric(paste(CCAT_chr16$areas)))

BCP_areas<-quantification(BCP_chr16,bedgraph_chr16)
BCP_chr16$areas<-paste(BCP_areas)
BCP_chr16$log_areas<-log10(as.numeric(paste(BCP_chr16$areas)))


SICER_areas<-quantification(SICER_chr16,bedgraph_chr16)
SICER_chr16$areas<-paste(SICER_areas)
SICER_chr16$log_areas<-log10(as.numeric(paste(SICER_chr16$areas)))


ZCL_areas<-quantification(ZCL_chr16,bedgraph_chr16)
ZCL_chr16$areas<-paste(ZCL_areas)
ZCL_chr16$log_areas<-log10(as.numeric(paste(ZCL_chr16$areas)))


areas_EXP1<-get(load("/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/areas_EXP1.Rdata"))
areas_EXP2<-get(load("/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/areas_EXP2.Rdata"))
areas_EXP3<-get(load("/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/06_SEÑALES_CHIPSEQ_ZCL_h3k4me3_ENE18/areas_exp3.Rdata"))

exp1_quantification_log<-log10(areas_EXP1)
exp2_quantification_log<-log10(areas_EXP2)
exp3_quantification_log<-log10(areas_EXP3)

plotter2<-data.frame('area' = NULL, 'Method' = NULL)
plotter2 <- rbind(plotter2, data.frame('area' = MACS_chr16$log_areas, 'Method' = rep('MACS2', nrow(MACS_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = JAMM_chr16$log_areas, 'Method' = rep('JAMM', nrow(JAMM_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = BAYES_chr16$log_areas, 'Method' = rep('BayesPeak', nrow(BAYES_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = CCAT_chr16$log_areas, 'Method' = rep('CCAT', nrow(CCAT_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = BCP_chr16$log_areas, 'Method' = rep('BCP', nrow(BCP_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = SICER_chr16$log_areas, 'Method' = rep('SICER', nrow(SICER_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = ZCL_chr16$log_areas, 'Method' = rep('PeakZCL', nrow(ZCL_chr16))))



pdf("violin_plot_areas_METODOS.pdf")
ggplot(plotter2, aes(x=Method, y=area)) + geom_violin() +geom_boxplot(width=0.05,outlier.shape = NA) +  theme_bw()
dev.off()

plotter2 <- rbind(plotter2, data.frame('area' = exp1_quantification_log, 'Method' = rep('Exp1', length(exp1_quantification_log))))
plotter2 <- rbind(plotter2, data.frame('area' = exp2_quantification_log, 'Method' = rep('Exp2', length(exp2_quantification_log))))
plotter2 <- rbind(plotter2, data.frame('area' = exp3_quantification_log, 'Method' = rep('Exp3', length(exp3_quantification_log))))

pdf("violin_plot_areas_ALL.pdf")
ggplot(plotter2, aes(x=Method, y=area)) + geom_violin() +geom_boxplot(width=0.05,outlier.shape = NA) +  theme_bw()
dev.off()


###########################################################
#########################################################

#EGR1

MACS_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_EGR1/MACS2/MACS2_EGR1_peaks.narrowPeak"
JAMM_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_EGR1/JAMM/selected_JAMM_EGR1.narrowPeak"
BAYES_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_EGR1/BAYES/bayes.bed"
RANGER_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_EGR1/RANGER/RANGER_EGR1.bed_region.bed"
SICER_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_EGR1/SICER/EGR1_chip_T-W200-G200-FDR0.01-island.bed"
ZCL_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_EGR1/ZCL/ZCL_peaks_g8_a1.2_k0_modified.bed"


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

require(ggplot2)

plotter <- data.frame('amplitude' = NULL, 'Method' = NULL)

plotter <- rbind(plotter, data.frame('amplitude' = MACS_chr16$amplitud, 'Method' = rep('MACS2', nrow(MACS_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = JAMM_chr16$amplitud, 'Method' = rep('JAMM', nrow(JAMM_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = BAYES_chr16$amplitud, 'Method' = rep('BayesPeak', nrow(BAYES_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = RANGER_chr16$amplitud, 'Method' = rep('RANGER', nrow(RANGER_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = SICER_chr16$amplitud, 'Method' = rep('SICER', nrow(SICER_chr16))))
plotter <- rbind(plotter, data.frame('amplitude' = ZCL_chr16$amplitud, 'Method' = rep('PeakZCL', nrow(ZCL_chr16))))

pdf("violin_plot_amplitudes_metodos_EGR1.pdf")
ggplot(plotter, aes(x=Method, y=amplitude)) + geom_violin() +geom_boxplot(width=0.02, outlier.shape = NA) +  theme_bw()
dev.off()

plotter <- rbind(plotter, data.frame('amplitude' = exp1$amplitud, 'Method' = rep('EXP1', nrow(exp1))))
plotter <- rbind(plotter, data.frame('amplitude' = exp2$amplitud, 'Method' = rep('EXP2', nrow(exp2))))
plotter <- rbind(plotter, data.frame('amplitude' = exp3$amplitud, 'Method' = rep('EXP3', nrow(exp3))))

pdf("violin_plot_amplitudes_ALL_EGR1.pdf")
ggplot(plotter, aes(x=Method, y=amplitude)) + geom_violin() +geom_boxplot(width=0.1, outlier.shape = NA) +  theme_bw()
dev.off()


bedgraph_chr16_file<-"/home/jgonzalez.69/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/cromosoma16_EGR1/chr16_genomecov.bedgraph"
bedgraph_chr16<-fread(bedgraph_chr16_file)
bedgraph_chr16<-as.data.frame(bedgraph_chr16)


MACS2_areas<-quantification(MACS_chr16,bedgraph_chr16)
MACS_chr16$areas<-paste(MACS2_areas)
MACS_chr16$log_areas<-log10(as.numeric(paste(MACS_chr16$areas)))

JAMM_areas<-quantification(JAMM_chr16,bedgraph_chr16)
JAMM_chr16$areas<-paste(JAMM_areas)
JAMM_chr16$log_areas<-log10(as.numeric(paste(JAMM_chr16$areas)))

BAYES_areas<-quantification(BAYES_chr16,bedgraph_chr16)
BAYES_chr16$areas<-paste(BAYES_areas)
BAYES_chr16$log_areas<-log10(as.numeric(paste(BAYES_chr16$areas)))

RANGER_areas<-quantification(RANGER_chr16,bedgraph_chr16)
RANGER_chr16$areas<-paste(RANGER_areas)
RANGER_chr16$log_areas<-log10(as.numeric(paste(RANGER_chr16$areas)))

SICER_areas<-quantification(SICER_chr16,bedgraph_chr16)
SICER_chr16$areas<-paste(SICER_areas)
SICER_chr16$log_areas<-log10(as.numeric(paste(SICER_chr16$areas)))

ZCL_areas<-quantification(ZCL_chr16,bedgraph_chr16)
ZCL_chr16$areas<-paste(ZCL_areas)
ZCL_chr16$log_areas<-log10(as.numeric(paste(ZCL_chr16$areas)))

plotter2<-data.frame('area' = NULL, 'Method' = NULL)
plotter2 <- rbind(plotter2, data.frame('area' = MACS_chr16$log_areas, 'Method' = rep('MACS2', nrow(MACS_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = JAMM_chr16$log_areas, 'Method' = rep('JAMM', nrow(JAMM_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = BAYES_chr16$log_areas, 'Method' = rep('BayesPeak', nrow(BAYES_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = RANGER_chr16$log_areas, 'Method' = rep('RANGER', nrow(RANGER_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = SICER_chr16$log_areas, 'Method' = rep('SICER', nrow(SICER_chr16))))
plotter2 <- rbind(plotter2, data.frame('area' = ZCL_chr16$log_areas, 'Method' = rep('PeakZCL', nrow(ZCL_chr16))))



pdf("violin_plot_areas_METODOS_EGR1.pdf")
ggplot(plotter2, aes(x=Method, y=area)) + geom_violin() +geom_boxplot(width=0.05,outlier.shape = NA) +  theme_bw()
dev.off()
