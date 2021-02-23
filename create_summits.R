
library(data.table)
library(compiler)
library(future)

enableJIT(3)
plan(multicore, workers=10)

options(future.globals.maxSize= +Inf)


señales<-list.files("/home/nostromo/data/pepe/BEDS_METODOS_EGR1/cromosomas_EGR1/", pattern=".chrom")
señales<-paste0("/home/nostromo/data/pepe/BEDS_METODOS_EGR1/cromosomas_EGR1/",señales)

señales_EGR1<-list.files("~/data/pepe/35_CHIPSEQ_EGR1_BAMS/chip/TMP/chip/", pattern=".chrom")
señales_EGR1<-paste0("~/data/pepe/35_CHIPSEQ_EGR1_BAMS/chip/TMP/chip/", señales_EGR1)

señales_h3k4me3<-list.files("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/cromosomas_h3k4me3/", pattern=".chrom")
señales_h3k4me3<-paste0("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/cromosomas_h3k4me3/",señales_h3k4me3)

señales_JUND<-list.files("/home/nostromo/data/pepe/37_JUND/chip/TMP/chip/", pattern=".chrom")
señales_JUND<-paste0("/home/nostromo/data/pepe/37_JUND/chip/TMP/chip/",señales_JUND)

señales_EGR1<-list.files("/home/nostromo/data/pepe/35_CHIPSEQ_EGR1_BAMS/chip/TMP/chip/", pattern=".chrom")
señales_EGR1<-paste0("/home/nostromo/data/pepe/35_CHIPSEQ_EGR1_BAMS/chip/TMP/chip/", señales_EGR1)

señales_SP1<-list.files("/home/nostromo/data/pepe/38_SP1/chip/TMP/chip/", pattern=".chrom")
señales_SP1<-paste0("/home/nostromo/data/pepe/38_SP1/chip/TMP/chip/",señales_SP1)


señales_h3k4me3<-list.files("/home/nostromo/data/pepe/36_h3k4me3/chip/TMP/chip/", pattern=".chrom")
señales_h3k4me3<-paste0("/home/nostromo/data/pepe/36_h3k4me3/chip/TMP/chip/",señales_h3k4me3)


cromosomas<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10", "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY", "chrM")

find_summits<-function(lista_señales,bed,lista_cromosomas){
    total_output<-data.frame()
    for (cromosoma in lista_cromosomas){
        cat("\n",cromosoma,"\n")
        temporal<-bed[which(bed$V1==cromosoma),]
        if(nrow(temporal)!=0){
        #output_chr1<-data.frame()
            temporal_output<-data.frame()
            lista_señales<-unlist(strsplit(lista_señales,".chrom"))
            #selected_signal<-lista_señales[which(strsplit(lista_señales,"_")==cromosoma)]
            selected_signal<-lista_señales[which(basename(lista_señales)==cromosoma)]
            selected_signal<-paste0(selected_signal,".chrom")
            cat("loading signal", selected_signal,"\n")
            selected_signal<-fread(selected_signal)
            selected_signal<-as.data.frame(selected_signal)
            for (i in 1:nrow(temporal)){
                cat(i,",")
                    if(temporal$V2[i]==0){
                        temporal$V2[i]<-1
                    }
                    tmp<-selected_signal[which(selected_signal$V2==temporal$V2[i]):which(selected_signal$V2==temporal$V3[i]),]
                    max<-tmp$V2[which(tmp$V3==max(tmp$V3))]
                    max<-max[1]
                    end_max<-max+1
                    name_summit<-paste0(temporal$V2[i],"_",temporal$V3[i])
                    out<-data.frame("chr"=cromosoma,"start"=max, "end"=end_max, "name"=name_summit)
                    temporal_output<-rbind(temporal_output,out)
                    

            }
            if(nrow(temporal_output)>0){
            total_output<-rbind(total_output,temporal_output)
            }else{
                next
            }
        }
    }
    return(total_output)
}

find_summits_c<-cmpfun(find_summits)

##### BAYES

Bayes<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_EGR1/BAYES/bayes.bed")
Bayes_h3k4me3<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/BAYES/BayesPeak.bed")


Bayes_summits<-find_summits_c(señales,Bayes, cromosomas)
write.table(Bayes_summits,"/home/nostromo/data/pepe/BEDS_METODOS_EGR1/BAYES/BayesPeak_EGR1_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

Bayes_summits_h3k4me3<-find_summits_c(señales_h3k4me3,Bayes_h3k4me3, cromosomas)
write.table(Bayes_summits_h3k4me3,"/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/BAYES/BayesPeak_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

Bayes_SP1<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/BAYESPEAK_RESULTS/BayesPeak_output_SP1.bed")
Bayes_summits_SP1<-find_summits_c(señales_SP1,Bayes_SP1, cromosomas)
write.table(Bayes_summits_SP1,"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/BAYESPEAK_RESULTS/BayesPeak_SP1_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

##### JAMM

JAMM<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_EGR1/JAMM/selected_JAMM_EGR1.narrowPeak")
JAMM_h3k4me3<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/JAMM/selected_JAMM.narrowPeak")

Jamm_summits_EGR1<-find_summits_c(señales, JAMM,cromosomas)

write.table(Jamm_summits_EGR1,"/home/nostromo/data/pepe/BEDS_METODOS_EGR1/JAMM/JAMM_EGR1_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

JAMM_summits_h3k4me3<-find_summits_c(señales_h3k4me3, JAMM_h3k4me3, cromosomas)
write.table(JAMM_summits_hek4me3,"/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/JAMM/JAMM_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")


##### SICER

SICER_EGR1<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_EGR1/SICER/EGR1_chip_T-W200-G200-FDR0.01-island.bed")
SICER_summits_EGR1<-find_summits_c(señales, SICER_EGR1,cromosomas)
write.table(SICER_summits_EGR1,"/home/nostromo/data/pepe/BEDS_METODOS_EGR1/SICER/SICER_EGR1_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")






SICER_h3k4me3<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/SICER/CHIP_h3k4me3-W200-G200-FDR0.01-island.bed")
SICER_summits_h3k4me3<-find_summits_c(señales_h3k4me3, SICER_h3k4me3,cromosomas)

write.table(SICER_summits_h3k4me3,"/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/SICER/SICER_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")


##### MACS2

MACS2_EGR1<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_EGR1/MACS2/MACS2_EGR1_peaks.narrowPeak")
MACS2_summits_EGR1<-find_summits_c(señales_SP1, MACS2_EGR1,cromosomas)
write.table(MACS2_summits_EGR1,"/home/nostromo/data/pepe/BEDS_METODOS_EGR1/MACS2/MACS2_EGR1_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")


#summit de MACS2 h3k4me3 está sacado del software

##### RANGER

RANGER_EGR1<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_EGR1/RANGER/RANGER_EGR1.bed_region.bed")
RANGER_EGR1<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/RANGER_RESULTS/RANGER_EGR1.bed_region.bed")

RANGER_SP1<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_SP1/RANGER/RANGER_SP1.bed_region.bed")
RANGER_SP1<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/RANGER_RESULTS/RANGER_SP1.bed_region.bed")
RANGER_JUND<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/RANGER_RESULTS/RANGER_JUND.bed_region.bed")
RANGER_JUND<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/RANGER_RESULTS/RANGER_JUND.bed_region.bed")
RANGER_h3k4me3<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/RANGER_RESULTS/RANGER_h3k4me3.bed_region.bed")



RANGER_summits_EGR1<-find_summits_c(señales_EGR1, RANGER_EGR1,cromosomas)
#write.table(RANGER_summits_EGR1,"/home/nostromo/data/pepe/BEDS_METODOS_EGR1/RANGER/RANGER_EGR1_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")
write.table(RANGER_summits_EGR1,"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/RANGER_RESULTS/RANGER_EGR1_summits_pepe.bed", row.names=F, col.names=F, quote=F, sep="\t")

RANGER_summits_SP1<-find_summits_c(señales_SP1, RANGER_SP1,cromosomas)
#write.table(RANGER_summits_SP1,"/home/nostromo/data/pepe/BEDS_METODOS_SP1/RANGER/RANGER_SP1_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")
write.table(RANGER_summits_SP1,"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/RANGER_RESULTS/RANGER_SP1_summits_pepe.bed", row.names=F, col.names=F, quote=F, sep="\t")

RANGER_summits_h3k4me3<-find_summits_c(señales_h3k4me3, RANGER_h3k4me3,cromosomas)
#write.table(RANGER_summits_h3k4me3,"/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/RANGER/RANGER_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")
write.table(RANGER_summits_h3k4me3,"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/RANGER_RESULTS/RANGER_h3k4me3_summits_pepe.bed", row.names=F, col.names=F, quote=F, sep="\t")

RANGER_summits_JUND<-find_summits_c(señales_JUND, RANGER_JUND,cromosomas)
#write.table(RANGER_summits_JUND,"/home/nostromo/data/pepe/BEDS_METODOS_JUND/RANGER/RANGER_JUND_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")
write.table(RANGER_summits_JUND,"/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/RANGER_RESULTS/RANGER_JUND_summits_pepe.bed", row.names=F, col.names=F, quote=F, sep="\t")




##### ZCL
ZCL_SP1<-read.table("/home/nostromo/data/pepe/38_SP1/chip/OUTPUT23/PeakZCL_peaks_zCros_24_clus_1_norm_SES_min_120_max_1000.sorted_final.bed")
ZCL_summits_SP1<-find_summits_c(señales_SP1, ZCL_SP1,cromosomas)
write.table(ZCL_summits_SP1,"/home/nostromo/data/pepe/38_SP1/chip/OUTPUT23/PeakZCL_peaks_zCros_24_clus_1_norm_SES_min_120_max_1000.sorted_final_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

ZCL_JUND<-read.table("/home/nostromo/data/pepe/37_JUND/chip/OUTPUT24/PeakZCL_peaks_zCros_24_clus_1_norm_SES_min_60_max_1000.sorted_final.bed")
ZCL_summits_JUND<-find_summits_c(señales_JUND, ZCL_JUND,cromosomas)
write.table(ZCL_summits_JUND,"/home/nostromo/data/pepe/37_JUND/chip/OUTPUT24/PeakZCL_peaks_zCros_24_clus_1_norm_SES_min_60_max_1000.sorted_final_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")



########################################3

ZCL_EGR1<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_EGR1/ZCL/ZCL_peaks.final.bed")
ZCL_summits_EGR1<-find_summits_c(señales, ZCL_EGR1,cromosomas)
write.table(ZCL_summits_EGR1,"/home/nostromo/data/pepe/BEDS_METODOS_EGR1/ZCL/ZCL_EGR1_summits_nueva.bed", row.names=F, col.names=F, quote=F, sep="\t")
###################33

ZCL_EGR1<-read.table("/home/nostromo/data/pepe/35_CHIPSEQ_EGR1_BAMS/chip/OUTPUT24/PeakZCL_peaks_zCros_22_clus_1_norm_SES_min_100_max_1500.sorted_final.bed")
ZCL_summits_EGR1<-find_summits_c(señales_EGR1, ZCL_EGR1,cromosomas)
write.table(ZCL_summits_EGR1,"/home/nostromo/data/pepe/35_CHIPSEQ_EGR1_BAMS/chip/OUTPUT24/ZCL_EGR1_summits_nueva.bed", row.names=F, col.names=F, quote=F, sep="\t")

ZCL_h3k4me3<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/ZCL/ZCL_peaks.bed")
ZCL_summits_h3k4me3<-find_summits_c(señales_h3k4me3, ZCL_h3k4me3,cromosomas)
write.table(ZCL_summits_h3k4me3,"/home/nostromo/data/pepe/36_h3k4me3/chip/OUTPUT23/ZCL_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

ZCL_h3k4me3<-read.table("/home/nostromo/data/pepe/36_h3k4me3/chip/OUTPUT23/PeakZCL_peaks_zCros_22_clus_10_norm_SES_min_150_max_6500.sorted_final.bed")
ZCL_summits_h3k4me3<-find_summits_c(señales_h3k4me3, ZCL_h3k4me3,cromosomas)
write.table(ZCL_summits_h3k4me3,"/home/nostromo/data/pepe/36_h3k4me3/chip/OUTPUT23/ZCL_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")


#### BCP

BCP_h3k4me3<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/BCP/h3k4me3_bcp.bed_region.bed")
BCP_summits_h3k4me3<-find_summits_c(señales_h3k4me3, BCP_h3k4me3,cromosomas)
write.table(BCP_summits_h3k4me3,"/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/BCP/BCP_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")


## BCP solo lo usamos para histonas
CCAT_h3k4me3<-read.table("~/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/CCAT/CCAT_h3k4me3.bed_region_passed.bed")
CCAT_summits_h3k4me3<-find_summits_c(señales_h3k4me3, CCAT_h3k4me3,cromosomas)
write.table(CCAT_summits_h3k4me3,"~/dato-activo/03_Analysis/jgonzalez69/ELI_ZCL/BEDS_METODOS_h3k4me3/CCAT/CCAT_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

### CCAT summits es de la herramienta y no lo usamos para EGR1 (no funcionaba el de la herramienta para lo q usa Eli y lo creo con mi script)

##### EXP1

EXP1<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/regiones_PEPE_h3k4me3_buena.bed", header=T)
EXP1_summits_h3k4me3<-find_summits_c(señales_h3k4me3, EXP1,cromosomas)

write.table(EXP1_summits_h3k4me3,"/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/EXP1_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

### EXP2

EXP2<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/regions_eli_todo.bed")
EXP2_summits_h3k4me3<-find_summits_c(señales_h3k4me3, EXP2,cromosomas)

write.table(EXP2_summits_h3k4me3,"/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/EXP2_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")

### EXP3

EXP3<-read.table("/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/CHIP_h3k4me3_T_regions_macarena.bed")
EXP3_summits_h3k4me3<-find_summits_c(señales_h3k4me3, EXP3,cromosomas)

write.table(EXP3_summits_h3k4me3,"/home/nostromo/data/pepe/BEDS_METODOS_h3k4me3/EXP3_h3k4me3_summits.bed", row.names=F, col.names=F, quote=F, sep="\t")
