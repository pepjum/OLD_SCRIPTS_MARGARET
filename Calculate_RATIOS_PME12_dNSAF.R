path<-"/home/margaret/data/pepe/10_PME12_SPECTRAL_COUNTING_8_5_19/Results_dNSAF/"

files<-list.files(path,pattern="dNSAF.txt")

files_w_path<-paste0(path,files)
options(stringsAsFactors = FALSE)
options(scipen=999)
for (i in 1: length(files_w_path)){
        cat(i,"\n")
        if (i ==1){
            name<-basename(files_w_path[i])
            name_splitted_tmp<-strsplit(name,"_dNSAF")[[1]][1]
            #name_splitted_2<-strsplit(name_splitted_tmp,"_dataPSM")[[1]][1]
            data_all<-read.csv2(files_w_path[i], header=T, sep="\t")
            data_all$sample<-name_splitted_tmp
        }else{
            name<-basename(files_w_path[i])
            name_splitted_tmp<-strsplit(name,"_dNSAF")[[1]][1]
            #name_splitted_2<-strsplit(name_splitted_tmp,"_dataPSM")[[1]][1]
            tmp<-read.csv2(files_w_path[i], header=T, sep="\t")
            tmp$sample<-name_splitted_tmp
            data_all<-rbind(data_all,tmp)
        }
}

names(data_all)[1]<-"ProteinAccession"
#eliminar decoys
#data_all<-data_all[!grepl("DECOY",data_all$ProteinAccession),]
#data_all$dNSAF<-format(round(as.numeric(data_all$dNSAF), 8), nsmall = 8)

matriz_all<-matrix(NA, length(unique(data_all$ProteinAccession)), length(unique(data_all$sample)))
colnames(matriz_all) <-unique(data_all$sample)
rownames(matriz_all) <- paste(unique(data_all$ProteinAccession))

#CREAR MATRIZ PARA CONTEAR EN CUANTAS MUESTRAS SE ESTA CUANTIFICANDO CADA PROTEINA

for(protein in unique(data_all$ProteinAccession)){
    print(protein)
    tmp <- data_all[data_all$ProteinAccession == protein, ]
    for (sample in unique(tmp$sample)){
        matriz_all[protein, sample] <- tmp[tmp$sample == sample, 'dNSAF']

    }
}

matriz_B<-matrix(NA, length(unique(data_all$ProteinAccession)), 2)
matriz_all<-cbind(matriz_all,matriz_B)
colnames(matriz_all)<-c(unique(data_all$sample),"n_of_samples_quantified_A","n_of_samples_quantified_B")

for(protein in unique(data_all$ProteinAccession)){
    print(protein)
    matriz_all[protein,"n_of_samples_quantified_A"]<-sum(!is.na(matriz_all[protein,][1:4]))
    matriz_all[protein,"n_of_samples_quantified_B"]<-sum(!is.na(matriz_all[protein,][5:8]))
}

df_all<-as.data.frame(matriz_all)
df_all$ProteinAccession<-rownames(matriz_all)

#dat<-df_all_filtered_by_nsamples_quantified
df_all$A_1<-as.numeric(df_all$A_1)
df_all$A_2<-as.numeric(df_all$A_2)
df_all$A_3<-as.numeric(df_all$A_3)
df_all$A_4<-as.numeric(df_all$A_4)
df_all$B_1<-as.numeric(df_all$B_1)
df_all$B_2<-as.numeric(df_all$B_2)
df_all$B_3<-as.numeric(df_all$B_3)
df_all$B_4<-as.numeric(df_all$B_4)

#df_all<-transform(df_all, media_dNSAF_A = rowMeans(df_all[1:4]), na.rm = TRUE)
#df_all$media_dNSAF_A<-apply(df_all[,1:4],1, function(x) mean(x))
df_all$media_dNSAF_A<- rowMeans(df_all[,1:4], na.rm=TRUE)
df_all$media_dNSAF_B<-rowMeans(df_all[,5:8], na.rm=TRUE)





for_output<-df_all
for_output$AC<-sapply(strsplit(paste(for_output$ProteinAccession),"\\|"),"[",2)

df_all$RATIO_RAW<-(df_all$media_dNSAF_B/df_all$media_dNSAF_A)

df_all$A_1_log<-log10(df_all$A_1)
df_all$A_2_log<-log10(df_all$A_2)
df_all$A_3_log<-log10(df_all$A_3)
df_all$A_4_log<-log10(df_all$A_4)
df_all$B_1_log<-log10(df_all$B_1)
df_all$B_2_log<-log10(df_all$B_2)
df_all$B_3_log<-log10(df_all$B_3)
df_all$B_4_log<-log10(df_all$B_4)

df_all$media_dNSAF_A_log<- rowMeans(df_all[,15:18], na.rm=TRUE)
df_all$media_dNSAF_B_log<-rowMeans(df_all[,19:22], na.rm=TRUE)

df_all$Ratio_log<- df_all$media_dNSAF_B_log - df_all$media_dNSAF_A_log
df_all_filtered_by_nsamples_quantified<-df_all[which((as.numeric(df_all$n_of_samples_quantified_A) >=1) & (as.numeric(df_all$n_of_samples_quantified_B) >=1)  ),]
dat<-(df_all_filtered_by_nsamples_quantified)
save(dat,file="quantification_Ratios_dNSAF_at_PSM_level_filtered_by_nsamples_quantified.Rdata")

dat[is.nan(dat$Ratio_log)]<-0

#### cargar los objetos de protFDR y quedarnos solo con las proteinas que han pasado todos los filtros)

path<-"/home/margaret/data/pepe/10_PME12_SPECTRAL_COUNTING_8_5_19/Dat_Files/"

files<-list.files(path,pattern="protFDR_Filter.rda")

files_w_path<-paste0(path,files)

for (i in 1: length(files_w_path)){
        cat(i,"\n")
        if (i ==1){
            results_PROTFDR<-get(load(files_w_path[i]))
        }else{
            tmp<-get(load(files_w_path[i]))
            results_PROTFDR<-rbind(results_PROTFDR,tmp)
        }
}

#delete DECOY results

dataPSMMat_filterAA_psmFDR_Filter <- results_PROTFDR[results_PROTFDR$database == "T",]

ProteinACC_of_proteins_ProtFDR_filtered<-unique(paste(dataPSMMat_filterAA_psmFDR_Filter$ProteinAccession))


#### filter proteins that not passed PROTFDR

object_final<-dat[(dat$ProteinAccession %in% ProteinACC_of_proteins_ProtFDR_filtered), ]

save(object_final, file="quantification_ratios_dNSAF_At_ProtFDR_level.Rdata")

#################################################################################################
# plots histograms of intensities separated by organism & % of wrongly quantified

object_final$organism<-lapply(strsplit(paste(object_final$ProteinAccession),"_"),"[",2)

HUMAN<-object_final[which(object_final$organism =="HUMAN"),]   #3597
ECO<-object_final[ grepl("^ECO",object_final$organism),] #632
YEAST<-object_final[grepl("^YE",object_final$organism),] #1417
YEASTV<-object_final[grepl("^S",object_final$organism),]

YEAST_TOT<-rbind(YEAST,YEASTV)

proteins_good_quantified_HUMAN<-HUMAN[which((HUMAN$Ratio_log >= -0.15) & (HUMAN$Ratio_log <=0.11)),] #1665

proteins_good_quantified_YEAST<-YEAST_TOT[which((YEAST_TOT$Ratio_log >=0.14) & (YEAST_TOT$Ratio_log <=0.42)),] #797

proteins_good_quantified_ECO<-ECO[which((ECO$Ratio_log >= -0.757) &(ECO$Ratio_log <=-0.488 )),] #19


####################### HUMAN


SPC_proteins_HUMAN<-proteins_good_quantified_HUMAN

#Dispersion of the ratios HUMAN
sd_HUMAN<-sd(SPC_proteins_HUMAN$RATIO_RAW)
mean_HUMAN<-mean(SPC_proteins_HUMAN$RATIO_RAW)

HUMAN_disp<-(sd_HUMAN/mean_HUMAN)*100

plot(SPC_proteins_HUMAN$RATIO_log10)


################# YEAST

SPC_proteins_YEAST<-proteins_good_quantified_YEAST

#dispersion
sd_YEAST<-sd(SPC_proteins_YEAST$RATIO_RAW)
mean_YEAST<-mean(SPC_proteins_YEAST$RATIO_RAW)

YEAST_disp<-(sd_YEAST/mean_YEAST)*100


################# ECOLI

SPC_proteins_ECO<-proteins_good_quantified_ECO

#dispersion
sd_ECO<-sd(SPC_proteins_ECO$RATIO_RAW)
mean_ECO<-mean(SPC_proteins_ECO$RATIO_RAW)

ECO_disp<-(sd_ECO/mean_ECO)*100

#plots

pdf("Distr_ratios_by_species_dNSAF.pdf")
par(mfrow=c(2,2))
hist(HUMAN$RATIO_RAW)
abline(v=c(0.7,1,1.3), col=c("red","blue","red"), lwd=2, lty=1)
hist(YEAST_TOT$RATIO_RAW)
abline(v=c(1.4,2,2.6), col=c("red","blue","red"), lwd=2, lty=1)
hist(ECO$RATIO_RAW)
abline(v=c(0.175,0.25,0.325), col=c("red","blue","red"), lwd=2, lty=1)
dev.off()

pdf("Distr_ratios_log10_by_species_dNSAF.pdf")
par(mfrow=c(2,2))
hist(HUMAN$Ratio_log)
abline(v=c(-0.15,0,0.11), col=c("red","blue","red"), lwd=2, lty=1)
hist(YEAST_TOT$Ratio_log)
abline(v=c(0.146,0.3,0.414), col=c("red","blue","red"), lwd=2, lty=1)
hist(ECO$Ratio_log)
abline(v=c(-0.757,-0.60,-0.488), col=c("red","blue","red"), lwd=2, lty=1)
dev.off()


tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio_log"=HUMAN$Ratio_log)
tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio_log"=ECO$Ratio_log)
tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST_TOT)),"ratio_log"=YEAST_TOT$Ratio_log)

library(ggplot2)

plotter<-data.frame("organism"=NULL,"ratio_log"=NULL)
plotter<-rbind(plotter,tmp_HUMAN)
plotter<-rbind(plotter,tmp_ECO)
plotter<-rbind(plotter,tmp_YEAST)

#scatter `plot`
pdf("box_plot_distr_ratios_log_by_species_dNSAF.pdf")
ggplot(data=plotter, aes(x=organism, y=ratio_log, fill=organism)) + geom_boxplot() + geom_hline(yintercept=c(-0.15,0,0.11,-0.757,-0.6,-0.488,0.146,0.3,0.414), color=c("red","red","red","green","green","green","blue","blue","blue"), linetype=c("dashed","solid","dashed","dashed","solid","dashed","dashed","solid","dashed"))
dev.off()

tmp_HUMAN_scatter<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio_log"=HUMAN$Ratio_log, "media_B"=HUMAN$media_dNSAF_B_log)
tmp_ECO_scatter<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio_log"=ECO$Ratio_log,"media_B"=ECO$media_dNSAF_B_log)
tmp_YEAST_scatter<-data.frame("organism"=rep("Yeast",nrow(YEAST_TOT)),"ratio_log"=YEAST_TOT$Ratio_log,"media_B"=YEAST_TOT$media_dNSAF_B_log)

plotter_scatter<-data.frame("organism"=NULL,"ratio_log"=NULL, "media_B"=NULL)
plotter_scatter<-rbind(plotter_scatter,tmp_HUMAN_scatter)
plotter_scatter<-rbind(plotter_scatter,tmp_ECO_scatter)
plotter_scatter<-rbind(plotter_scatter,tmp_YEAST_scatter)




pdf("scatter_plot_distr_ratios_log_by_species_dNSAF.pdf")
ggplot(data=plotter_scatter, aes(x=media_B, y=ratio_log, fill=organism, color=organism)) + geom_point() +  geom_hline(yintercept=c(-0.15,0,0.11,-0.757,-0.6,-0.488,0.146,0.3,0.414), color=c("red","red","red","green","green","green","blue","blue","blue"), linetype=c("dashed","solid","dashed","dashed","solid","dashed","dashed","solid","dashed"))
dev.off()


source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
library(ggplot2)
#PLOTS DETOUTLIERS

subset_ELI_no_N_dNSAF<-object_final[,c(15:22,11)]

acol = sample(brewer.pal(8, "Dark2"), ncol(subset_ELI_no_N_dNSAF)-1, replace = (8 < ncol(subset_ELI_no_N_dNSAF)-1))
xlim = c(min(na.omit(subset_ELI_no_N_dNSAF[,c(1:8)])),max(na.omit(subset_ELI_no_N_dNSAF[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(subset_ELI_no_N_dNSAF[,c(1:8)])),method="average")

pdf(file="QC_PME12_dNSAF_May19.pdf", width = 15, height = 15)
boxplot(subset_ELI_no_N_dNSAF[,c(1:8)],names=colnames(subset_ELI_no_N_dNSAF)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(subset_ELI_no_N_dNSAF[,c(1:8)]), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_no_N_dNSAF[,c(1:8)], tipo, "DetOut_PME12_dNSAF_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_no_N_dNSAF[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_no_N_dNSAF[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_dNSAF_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

#### epecies separadas

HUMAN<-object_final[which(object_final$organism =="HUMAN"),]   #3597
ECO<-object_final[ grepl("^ECO",object_final$organism),] #632
YEAST<-object_final[grepl("^YE",object_final$organism),] #1424
YEASTV<-object_final[grepl("^S",object_final$organism),]
YEAST<-rbind(YEAST,YEASTV)

HUMAN_ELI_no_N_dNSAF<-HUMAN[,c(15:22,11)]
ECO_ELI_no_N_dNSAF<-ECO[,c(15:22,11)]
YEAST_ELI_no_N_dNSAF<-YEAST[,c(15:22,11)]

####### HUMAN
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesVikv2.R")

acol = sample(brewer.pal(8, "Dark2"), ncol(HUMAN_ELI_no_N_dNSAF)-1, replace = (8 < ncol(HUMAN_ELI_no_N_dNSAF)-1))
xlim = c(min(na.omit(HUMAN_ELI_no_N_dNSAF[,c(1:8)])),max(na.omit(HUMAN_ELI_no_N_dNSAF[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(HUMAN_ELI_no_N_dNSAF[,c(1:8)])),method="average")

pdf(file="QC_PME12_dNSAF_HUMAN_May19.pdf", width = 15, height = 15)
boxplot(HUMAN_ELI_no_N_dNSAF[,c(1:8)],names=colnames(HUMAN_ELI_no_N_dNSAF)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(HUMAN_ELI_no_N_dNSAF[,c(1:8)]), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(HUMAN_ELI_no_N_dNSAF[,c(1:8)], tipo, "DetOut_PME12_dNSAF_HUMAN_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(HUMAN_ELI_no_N_dNSAF[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(HUMAN_ELI_no_N_dNSAF[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_dNSAF_HUMAN_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

#### ECOLI

acol = sample(brewer.pal(8, "Dark2"), ncol(ECO_ELI_no_N_dNSAF)-1, replace = (8 < ncol(ECO_ELI_no_N_dNSAF)-1))
xlim = c(min(na.omit(ECO_ELI_no_N_dNSAF[,c(1:8)])),max(na.omit(ECO_ELI_no_N_dNSAF[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(ECO_ELI_no_N_dNSAF[,c(1:8)])),method="average")

pdf(file="QC_PME12_dNSAF_ECO_May19.pdf", width = 15, height = 15)
boxplot(ECO_ELI_no_N_dNSAF[,c(1:8)],names=colnames(ECO_ELI_no_N_dNSAF)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(ECO_ELI_no_N_dNSAF[,c(1:8)]), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(ECO_ELI_no_N_dNSAF[,c(1:8)], tipo, "DetOut_PME12_dNSAF_ECO_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(ECO_ELI_no_N_dNSAF[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(ECO_ELI_no_N_dNSAF[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_dNSAF_ECO_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

### YEAST_

acol = sample(brewer.pal(8, "Dark2"), ncol(YEAST_ELI_no_N_dNSAF)-1, replace = (8 < ncol(YEAST_ELI_no_N_dNSAF)-1))
xlim = c(min(na.omit(YEAST_ELI_no_N_dNSAF[,c(1:8)])),max(na.omit(YEAST_ELI_no_N_dNSAF[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(YEAST_ELI_no_N_dNSAF[,c(1:8)])),method="average")

pdf(file="QC_PME12_dNSAF_YEAST_May19.pdf", width = 15, height = 15)
boxplot(YEAST_ELI_no_N_dNSAF[,c(1:8)],names=colnames(YEAST_ELI_no_N_dNSAF)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(YEAST_ELI_no_N_dNSAF[,c(1:8)]), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_YEAST)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(YEAST_ELI_no_N_dNSAF[,c(1:8)], tipo, "DetOut_PME12_dNSAF_YEAST_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(YEAST_ELI_no_N_dNSAF[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(YEAST_ELI_no_N_dNSAF[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_dNSAF_YEAST_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()





 path_files<-"/home/margaret/data/pepe/10_PME12_SPECTRAL_COUNTING_8_5_19/Dat_Files/"
 name_sample_files<-list.files(path=path_files, pattern="dataPSMMat_filterAA_PSMFDR_filter.rda")
 file_path<-paste0(path_files,name_sample_files)


for (i in 1: length(file_path)){
        cat(i,"\n")
        if (i ==1){
            data_all<-get(load(files_w_path[i]))
        }else{
            tmp<-get(load(files_w_path[i]))
            data_all<-rbind(data_all,tmp)
        }
}
