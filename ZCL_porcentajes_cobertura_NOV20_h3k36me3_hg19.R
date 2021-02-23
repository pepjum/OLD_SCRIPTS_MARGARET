



#referencia<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3/hg19/E001-H3K36me3_sorted.broadPeak") #estaria mal. cambia por la remapped
referencia<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/remapped_homo_sapiens.GRCh38.K562.H3K36me3.ccat_histone.peaks.20190329_hg19.bed")
chromsizes<-read.table("~/data/00_References/bowtieIndexes/hg19_Bowtie2/hg19.chrom.sizes")



product<-function(referencia, chrmosizes, metodo){
    totales<-c()
    totales_ref<-c()
    cromosomas<-c("chr1","chr2","chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17", "chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")
    for(chr in cromosomas){
        cat(chr,"\n")
        vector_referencia<-rep(0,chromsizes$V2[which(chromsizes$V1==chr)],)
        vector_metodo<-rep(0,chromsizes$V2[which(chromsizes$V1==chr)],)
        metodo_chr<-metodo[which(metodo$V1==chr),]
        referencia_chr<-referencia[which(referencia$V1==chr),]
        if(nrow(metodo_chr) >0){
            #cat(" vector metodo","\n")
            for(i in 1:nrow(metodo_chr)){
                region_start<-metodo_chr[i,2]
                region_end<-metodo_chr[i,3]
                vector_metodo[region_start:region_end]<-1
            }
            #cat(" vector referencia","\n")
            if(nrow(referencia_chr)>0){
                for(i in 1:nrow(referencia_chr)){
                    region_start<-referencia_chr[i,2]
                    region_end<-referencia_chr[i,3]
                    vector_referencia[region_start:region_end]<-1

                }
            }else{
                vector_referencia<-vector_referencia

            }
            sum_referencia_chr<-sum(vector_referencia)
            #cat ("producto",chr,"\n")
            total_chr<-sum(vector_referencia * vector_metodo)
            cat("solapamiento ", chr," ",total_chr,"\n")
            if (!is.na(total_chr)){
            totales<-c(totales, total_chr)
            totales_ref<-c(totales_ref, sum_referencia_chr)
            cat(totales,"\n")
            cat(totales_ref,"\n")
            }
        }
    }
    newlist<-list("metodo"=totales,"referencia"= totales_ref)
    return(newlist)
}

#### MACS

#MACS<-read.table("/home/nostromo/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3_hg19/MACS2/MACS2_h3k36me3_hg19_peaks.broadPeak")
MACS<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/MACS2/MACS2_h3k36me3_peaks.broadPeak")
MACS_solapamiento<-product(referencia, chromsizes ,MACS)
porcentaje_macs<-sum(MACS_solapamiento$metodo)/sum(MACS_solapamiento$referencia)*100

### SICER

SICER<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/SICER/CHIP_h3k36me3_hg19_NOV20_sorted-W200-G600-islands-summary-FDR0.1.bed")
SICER_solapamiento<-product(referencia, chromsizes ,SICER)
porcentaje_SICER<-sum(SICER_solapamiento$metodo)/sum(SICER_solapamiento$referencia)*100

### BCP
BCP<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/BCP_0.01/BCP_0.01_h3k36me3_hg19_NOV.bed_region.bed")
BCP_solapamiento<-product(referencia, chromsizes ,BCP)
porcentaje_BCP<-sum(BCP_solapamiento$metodo)/sum(BCP_solapamiento$referencia)*100

### PeakZCL
PEAKZCL<-read.table("/home/nostromo/data/pepe/39_h3k36me3_hg19/chip/OUTPUT26/PeakZCL_peaks_zCros_4_clus_10_norm_SES_min_100_max_10000000.sorted_final.bed")
PEAKZCL_solapamiento<-product(referencia, chromsizes ,PEAKZCL)
porcentaje_PEAKZCL<-sum(PEAKZCL_solapamiento$metodo)/sum(PEAKZCL_solapamiento$referencia)*100

#PEAKZCL<-read.table("/home/nostromo/data/pepe/39_h3k36me3_hg19/chip/OUTPUT6/PeakZCL_peaks_zCros_6_clus_10_norm_SES_min_200_max_10000.sorted_final.bed")


#### CCAT

CCAT<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/CCAT/CCAT_h3k36me3_hg19_NOV20.bed_region.bed")
CCAT_solapamiento<-product(referencia, chromsizes ,CCAT)
porcentaje_CCAT<-sum(CCAT_solapamiento$metodo)/sum(CCAT_solapamiento$referencia)*100


### BayesPeak
BayesPeak<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/BayesPeak_h3k36me3_hg19_NOV20.bed")
BayesPeak_solapamiento<-product(referencia, chromsizes, BayesPeak)
porcentaje_BayesPeak<-sum(BayesPeak_solapamiento$metodo)/sum(BayesPeak_solapamiento$referencia)*100


## JAMM
JAMM<-read.table("/home/nostromo/data/pepe/19_CHIPSEQ_REVIEW_SIGNALS_PREPARATIONDIC19/h3k36me3_NOV20_hg19/JAMM/peaks/selected_JAMM_h3k36me3_NOV20_hg19.narrowPeak")

JAMM_solapamiento<-product(referencia, chromsizes, JAMM)
porcentaje_JAMM<-sum(JAMM_solapamiento$metodo)/sum(JAMM_solapamiento$referencia)*100

porcentajes_df<-data.frame("metodos"=c("MACS","SICER","BCP","CCAT","BayesPeak","JAMM","PEAKZCL"),"porcentaje"=c(porcentaje_macs,porcentaje_SICER,porcentaje_BCP,porcentaje_CCAT,porcentaje_BayesPeak,porcentaje_JAMM,porcentaje_PEAKZCL))
write.table(porcentajes_df, file="/home/nostromo/data/pepe/39_h3k36me3_hg19/porcentajes_cobertura_h3k36me3_hg19_DIC20.txt", col.names=T, sep="\t", quote=F)






############################# VUELVO a alinear h3k36me3 

#cat SRR40 SRR41 > CHIP_h3k36me3_NOV20.fastq

bowtie2 -p 2 -x ~/data/00_References/bowtieIndexes/hg19_Bowtie2/hg19 -U file1 -S CHIP_h3k36me3_hg19.sam


