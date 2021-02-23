
file<-read.csv2("/home/margaret/data/pepe/10_PME12_SPECTRAL_COUNTING_8_5_19/Results_MaxQuant/LFQ_24_5_2019_MaxQuant.txt", header=T, sep="\t")


File_arr<-file[-c(1,2),]
protein_grups<-read.csv2("/home/margaret/data/pepe/10_PME12_SPECTRAL_COUNTING_8_5_19/Results_MaxQuant/proteinGroups.txt", header=T, sep="\t")

#process MaxQuant output to see psm and proteins used for SC


b<-unlist(protein_grups$Peptide.IDs)

exit_list<-c()
for (i in 1:length(b)){
    print(i)
    char_b<-as.character(b[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    for( j in 1:length(char_splitted))
    exit_list<-c(exit_list, char_splitted[j])
}

number_of_peptides_to_SC<-length(unique(paste(exit_list)))

psms<-unlist(protein_grups$MS.MS.IDs)

exit_psm<-c()

for (i in 1:length(psms)){
    print(i)
    char_b<-as.character(psms[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    for( j in 1:length(char_splitted))
    exit_psm<-c(exit_psm, char_splitted[j])
}

number_of_psm_used<-length(unique(paste(exit_psm)))

#filtrar protein groups que tengan organismos mezclados

proteins_id<-paste(unlist(protein_grups$Protein.IDs))

only_1_specie_in_each_protein_group<-c()

for (i in 1:length(proteins_id)){
    print(i)
    char_b<-as.character(proteins_id[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    protein<-paste(char_splitted[1])
    specie_to_search<-sapply(strsplit(protein,"_"),"[", 2)
    start_specie<-substr(specie_to_search, 1, 3)
    species_to_merge<-c()
    for( j in 1:length(char_splitted)){
    specie_to_merge<-sapply(strsplit(char_splitted[j],"_"),"[", 2)
    species_to_merge<-c(species_to_merge,specie_to_merge)
    }
    specie_to_merge_start<-c()
    for(k in 1:length(species_to_merge)){
        specie_to_merge_start[k]<-substr(species_to_merge[k], 1, 3)
        if(start_specie==specie_to_merge_start[k]){

            only_1_specie_in_each_protein_group[i]<-"TRUE"

        }else{
            only_1_specie_in_each_protein_group[i]<-"FALSE"

        }
    }
}

#### quitar groous con más de 1 specie

protein_grups$only_one_specie<-paste(only_1_specie_in_each_protein_group)

protein_grups_sin<-protein_grups[-grepl("REV_", paste(protein_grups$Protein.IDs)),]
data_final_selected<-protein_grups[which(protein_grups$only_one_specie=="TRUE"),]

#### separar por organismo
library(stringr)

protein_first<-paste(data_final_selected$Protein.IDs)

organism<-c()
for (i in 1:length(protein_first)){
    count_commas<-str_count(protein_first[i],";")
    if(count_commas==0){
        proteinName<-sapply(strsplit(protein_first[i],"\\|"),"[",3)
        organism_detected<-strsplit(proteinName,"_")[[1]][2]
        organism<-c(organism,organism_detected)
    }else{
        accession_selected<-strsplit(protein_first[i],";")[[1]][1]
        proteinName<-sapply(strsplit(accession_selected,"\\|"),"[",3)
        organism_detected<-strsplit(proteinName,"_")[[1]][2]
        organism<-c(organism,organism_detected)

    }
}
data_final_selected$organism<-paste(organism)

HUMAN<-data_final_selected[which(data_final_selected$organism =="HUMAN"),]   #3781
ECO<-data_final_selected[ grepl("^ECO",data_final_selected$organism),] #2211
YEAST<-data_final_selected[grepl("^YE",data_final_selected$organism),]
YEASTV<-data_final_selected[grepl("^S",data_final_selected$organism),]
YEAST_TOT<-rbind(YEAST,YEASTV)

####### psm por especie
psm_human<-paste(HUMAN$MS.MS.IDs)
psm_eco<-paste(ECO$MS.MS.IDs)
psm_yeast<-paste(YEAST_TOT$MS.MS.IDs)

exit_psm_human<-c()

for (i in 1:length(psm_human)){
    print(i)
    char_b<-as.character(psm_human[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    for( j in 1:length(char_splitted))
    exit_psm_human<-c(exit_psm_human, char_splitted[j])
}

number_of_psm_used_human<-length(unique(paste(exit_psm_human)))

exit_psm_eco<-c()

for (i in 1:length(psm_eco)){
    print(i)
    char_b<-as.character(psm_eco[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    for( j in 1:length(char_splitted))
    exit_psm_eco<-c(exit_psm_eco, char_splitted[j])
}

number_of_psm_used_eco<-length(unique(paste(exit_psm_eco)))

exit_psm_yeast<-c()

for (i in 1:length(psm_yeast)){
    print(i)
    char_b<-as.character(psm_yeast[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    for( j in 1:length(char_splitted))
    exit_psm_yeast<-c(exit_psm_yeast, char_splitted[j])
}

number_of_psm_used_yeast<-length(unique(paste(exit_psm_yeast)))

#peptidos por especie

b_human<-unlist(HUMAN$Peptide.IDs)
b_eco<-unlist(ECO$Peptide.IDs)
b_yeast<-unlist(YEAST_TOT$Peptide.IDs)


exit_list_human<-c()
for (i in 1:length(b_human)){
    print(i)
    char_b<-as.character(b_human[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    for( j in 1:length(char_splitted))
    exit_list_human<-c(exit_list_human, char_splitted[j])
}

number_of_peptides_to_SC_human<-length(unique(paste(exit_list_human)))

exit_list_eco<-c()
for (i in 1:length(b_eco)){
    print(i)
    char_b<-as.character(b_eco[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    for( j in 1:length(char_splitted))
    exit_list_eco<-c(exit_list_eco, char_splitted[j])
}

number_of_peptides_to_SC_eco<-length(unique(paste(exit_list_eco)))

exit_list_yeast<-c()
for (i in 1:length(b_yeast)){
    print(i)
    char_b<-as.character(b_yeast[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    for( j in 1:length(char_splitted))
    exit_list_yeast<-c(exit_list_yeast, char_splitted[j])
}

number_of_peptides_to_SC_yeast<-length(unique(paste(exit_list_yeast)))

#################### processing Perseus output

protein_groups<-unlist(File_arr$Protein.IDs)

only_1_specie_in_each_protein_group<-c()

for (i in 1:length(protein_groups)){
    print(i)
    char_b<-as.character(protein_groups[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    protein<-paste(char_splitted[1])
    specie_to_search<-sapply(strsplit(protein,"_"),"[", 2)
    start_specie<-substr(specie_to_search, 1, 3)
    species_to_merge<-c()
    for( j in 1:length(char_splitted)){
    specie_to_merge<-sapply(strsplit(char_splitted[j],"_"),"[", 2)
    species_to_merge<-c(species_to_merge,specie_to_merge)
    }
    specie_to_merge_start<-c()
    for(k in 1:length(species_to_merge)){
        specie_to_merge_start[k]<-substr(species_to_merge[k], 1, 3)
        if(start_specie==specie_to_merge_start[k]){

            only_1_specie_in_each_protein_group[i]<-"TRUE"

        }else{
            only_1_specie_in_each_protein_group[i]<-"FALSE"

        }
    }
}

#### quitar groous con más de 1 specie

File_arr$only_one_specie<-paste(only_1_specie_in_each_protein_group)


data_final_selected<-File_arr[which(File_arr$only_one_specie=="TRUE"),]

specie_final<-c()
for (i in 1:nrow(data_final_selected)){
    print(i)
    char_b<-paste(as.character(data_final_selected$Protein.IDs[i]))
    char_splitted<-strsplit(char_b,";")
    protein<-char_splitted[[1]][1]
    specie_detected<-sapply(strsplit(protein,"_"),"[", 2)
    specie_final<-c(specie_final,specie_detected)
}
data_final_selected$ESPECIE<-paste(specie_final)

data_final_selected[,1]<-as.numeric(paste(data_final_selected[,1]))
data_final_selected[,2]<-as.numeric(paste(data_final_selected[,2]))
data_final_selected[,3]<-as.numeric(paste(data_final_selected[,3]))
data_final_selected[,4]<-as.numeric(paste(data_final_selected[,4]))
data_final_selected[,5]<-as.numeric(paste(data_final_selected[,5]))
data_final_selected[,6]<-as.numeric(paste(data_final_selected[,6]))
data_final_selected[,7]<-as.numeric(paste(data_final_selected[,7]))
data_final_selected[,8]<-as.numeric(paste(data_final_selected[,8]))


data_final_selected$media_LFQ_A<-rowMeans(data_final_selected[,1:4])
data_final_selected$media_LFQ_B<-rowMeans(data_final_selected[,5:8])

data_final_selected$RATIO_log2<-(data_final_selected$media_LFQ_B - data_final_selected$media_LFQ_A)


save(data_final_selected,file="quantification_Ratios_MAXQUANT_PERSEUS_LFQ_at_PSM_level.Rdata")


################ disgregar el objecto por especies

decoys<-data_final_selected[grepl("^REV", data_final_selected$Protein.IDs),]

HUMAN<-data_final_selected[which(data_final_selected$ESPECIE =="HUMAN"),]   #2391
ECO<-data_final_selected[ grepl("^ECO",data_final_selected$ESPECIE),] #131
YEAST<-data_final_selected[grepl("^YE",data_final_selected$ESPECIE),] #815
YEASTV<-data_final_selected[grepl("^S",data_final_selected$ESPECIE),]

YEAST<-rbind(YEAST,YEASTV)

# one_peptide_LFQ<-data_final_selected[which(data_final_selected$Peptides==1),]
# two_or_more_peptides_LFQ<-data_final_selected[which(paste(data_final_selected$Peptides)>=2),]
# HUMAN_two_peptide<-HUMAN[which(paste(HUMAN$Peptides)>=2),]
# ECO_two_peptide<-ECO[which(paste(ECO$Peptides)>=2),]
# YEAST_two_peptide<-YEAST[which(paste(YEAST$Peptides)>=2),]
# HUMAN_bck<-HUMAN
# YEAST_bck<-YEAST
# ECO_bck<-ECO
# HUMAN<-HUMAN_two_peptide
# YEAST<-YEAST_two_peptide
# ECO<-ECO_two_peptide

proteins_good_quantified_HUMAN<-HUMAN[which((HUMAN$RATIO_log2 >=-0.514) & (HUMAN$RATIO_log2 <=0.378)),] #2065

proteins_good_quantified_YEAST<-YEAST[which((YEAST$RATIO_log2 >=0.485) & (YEAST$RATIO_log2 <=1.378)),] #653

proteins_good_quantified_ECO<-ECO[which((ECO$RATIO_log2 >= -2.514) &(ECO$RATIO_log2 <=-1.621 )),] #78


##### desviacion estandar  y dispersion
HUMAN$RATIO_RAW<-2^HUMAN$RATIO_log2
ECO$RATIO_RAW<-2^ECO$RATIO_log2
YEAST$RATIO_RAW<-2^YEAST$RATIO_log2


sd_HUMAN<-sd(HUMAN$RATIO_RAW)
mean_HUMAN<-mean(HUMAN$RATIO_RAW)

HUMAN_disp<-(sd_HUMAN/mean_HUMAN)*100

sd_ECO<-sd(ECO$RATIO_RAW)
mean_ECO<-mean(ECO$RATIO_RAW)

ECO_disp<-(sd_ECO/mean_ECO)*100

sd_YEAST<-sd(YEAST$RATIO_RAW)
mean_YEAST<-mean(YEAST$RATIO_RAW)

YEAST_disp<-(sd_YEAST/mean_YEAST)*100

pdf("Distr_ratios_by_species_LFQ.pdf")
par(mfrow=c(2,2))
hist(HUMAN$RATIO_RAW)
abline(v=c(0.7,1,1.3), col=c("red","blue","red"), lwd=2, lty=1)
hist(YEAST$RATIO_RAW)
abline(v=c(1.4,2,2.6), col=c("red","blue","red"), lwd=2, lty=1)
hist(ECO$RATIO_RAW)
abline(v=c(0.175,0.25,0.325), col=c("red","blue","red"), lwd=2, lty=1)
dev.off()

pdf("Distr_ratios_log2_by_species_LFQ.pdf")
par(mfrow=c(2,2))
hist(HUMAN$RATIO_log2)
abline(v=c(-0.514,0,0.378), col=c("red","blue","red"), lwd=2, lty=1)
hist(YEAST$RATIO_log2)
abline(v=c(0.485,1,1.378), col=c("red","blue","red"), lwd=2, lty=1)
hist(ECO$RATIO_log2)
abline(v=c(-2.514,-2,-1.621), col=c("red","blue","red"), lwd=2, lty=1)
dev.off()

tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio_log"=HUMAN$RATIO_log2)
tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio_log"=ECO$RATIO_log2)
tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST)),"ratio_log"=YEAST$RATIO_log2)

library(ggplot2)

plotter<-data.frame("organism"=NULL,"ratio_log"=NULL)
plotter<-rbind(plotter,tmp_HUMAN)
plotter<-rbind(plotter,tmp_ECO)
plotter<-rbind(plotter,tmp_YEAST)

#box `plot`
pdf("box_plot_distr_ratios_log_by_species_LFQ.pdf")
ggplot(data=plotter, aes(x=organism, y=ratio_log, fill=organism)) + geom_boxplot() + geom_hline(yintercept=c(-2.514,-2,-1.621,-0.514,0,0.378,0.485,1,1.378), color=c("green","green","green","red","red","red","blue","blue","blue"), linetype=c("dashed","solid","dashed","dashed","solid","dashed","dashed","solid","dashed"))
dev.off()

# tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio"=HUMAN$RATIO_RAW)
# tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio"=ECO$RATIO_RAW)
# tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST)),"ratio"=YEAST$RATIO_RAW)
#
# library(ggplot2)
#
# plotter<-data.frame("organism"=NULL,"ratio"=NULL)
# plotter<-rbind(plotter,tmp_HUMAN)
# plotter<-rbind(plotter,tmp_ECO)
# plotter<-rbind(plotter,tmp_YEAST)
#
# #box `plot`
# pdf("box_plot_distr_ratios_by_species_LFQ.pdf")
# ggplot(data=plotter, aes(x=organism, y=ratio, fill=organism)) + geom_boxplot() + geom_hline(yintercept=c(0.25,1,2), color=c("green","red","blue"))
# dev.off()

#scatter_plot

tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio_log"=HUMAN$RATIO_log2, "media_LFQ_B"=HUMAN$media_LFQ_B)
tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio_log"=ECO$RATIO_log2, "media_LFQ_B"=ECO$media_LFQ_B)
tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST)),"ratio_log"=YEAST$RATIO_log2, "media_LFQ_B"=YEAST$media_LFQ_B)



plotter<-data.frame("organism"=NULL,"ratio_log"=NULL, "media_LFQ_B"=NULL)
plotter<-rbind(plotter,tmp_HUMAN)
plotter<-rbind(plotter,tmp_ECO)
plotter<-rbind(plotter,tmp_YEAST)



#scatter `plot`
pdf("scatter_plot_distr_ratios_log2_by_species_LFQ.pdf")
ggplot(data=plotter, aes(x=media_LFQ_B, y=ratio_log, fill=organism, color=organism)) + geom_point() + geom_hline(yintercept=c(-2.514,-2,-1.621,-0.514,0,0.378,0.485,1,1.378), color=c("green","green","green","red","red","red","blue","blue","blue"),
linetype=c("dashed","solid","dashed","dashed","solid","dashed","dashed","solid","dashed"))
dev.off()

################333

#PLOTS DETOUTLIERS
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesVikv2.R")
#from data_final_selected
subset_ELI_no_N_LFQ<-data_final_selected[,c(1:8,23)]
#subset_ELI_no_N_LFQ<-one_peptide_LFQ[,c(1:8,23)]
#subset_ELI_no_N_LFQ<-two_or_more_peptides[,c(1:8,23)]

acol = sample(brewer.pal(8, "Dark2"), ncol(subset_ELI_no_N_LFQ)-1, replace = (8 < ncol(subset_ELI_no_N_LFQ)-1))
xlim = c(min(na.omit(subset_ELI_no_N_LFQ[,c(1:8)])),max(na.omit(subset_ELI_no_N_LFQ[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(subset_ELI_no_N_LFQ[,c(1:8)])),method="average")

pdf(file="QC_PME12_LFQ_May19.pdf", width = 15, height = 15)
boxplot(subset_ELI_no_N_LFQ[,c(1:8)],names=colnames(subset_ELI_no_N_LFQ)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(subset_ELI_no_N_LFQ[,c(1:8)]+1), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_no_N_LFQ[,c(1:8)], tipo, "DetOut_PME12_LFQ_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_no_N_LFQ[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_no_N_LFQ[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_LFQ_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

#### epecies separadas

HUMAN_ELI_no_N_LFQ<-HUMAN[,c(1:8,23)]
ECO_ELI_no_N_LFQ<-ECO[,c(1:8,23)]
YEAST_ELI_no_N_LFQ<-YEAST[,c(1:8,23)]

####### HUMAN

acol = sample(brewer.pal(8, "Dark2"), ncol(HUMAN_ELI_no_N_LFQ)-1, replace = (8 < ncol(HUMAN_ELI_no_N_LFQ)-1))
xlim = c(min(na.omit(HUMAN_ELI_no_N_LFQ[,c(1:8)])),max(na.omit(HUMAN_ELI_no_N_LFQ[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(HUMAN_ELI_no_N_LFQ[,c(1:8)])),method="average")

pdf(file="QC_PME12_LFQ_HUMAN_May19.pdf", width = 15, height = 15)
boxplot(HUMAN_ELI_no_N_LFQ[,c(1:8)],names=colnames(HUMAN_ELI_no_N_LFQ)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(HUMAN_ELI_no_N_LFQ[,c(1:8)]+1), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(HUMAN_ELI_no_N_LFQ[,c(1:8)], tipo, "DetOut_PME12_LFQ_HUMAN_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(HUMAN_ELI_no_N_LFQ[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(HUMAN_ELI_no_N_LFQ[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_LFQ_HUMAN_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

#### ECOLI

acol = sample(brewer.pal(8, "Dark2"), ncol(ECO_ELI_no_N_LFQ)-1, replace = (8 < ncol(ECO_ELI_no_N_LFQ)-1))
xlim = c(min(na.omit(ECO_ELI_no_N_LFQ[,c(1:8)])),max(na.omit(ECO_ELI_no_N_LFQ[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(ECO_ELI_no_N_LFQ[,c(1:8)])),method="average")

pdf(file="QC_PME12_LFQ_ECO_May19.pdf", width = 15, height = 15)
boxplot(ECO_ELI_no_N_LFQ[,c(1:8)],names=colnames(ECO_ELI_no_N_LFQ)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(ECO_ELI_no_N_LFQ[,c(1:8)]+1), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(ECO_ELI_no_N_LFQ[,c(1:8)], tipo, "DetOut_PME12_LFQ_ECO_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(ECO_ELI_no_N_LFQ[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(ECO_ELI_no_N_LFQ[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_LFQ_ECO_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

### YEAST_

acol = sample(brewer.pal(8, "Dark2"), ncol(YEAST_ELI_no_N_LFQ)-1, replace = (8 < ncol(YEAST_ELI_no_N_LFQ)-1))
xlim = c(min(na.omit(YEAST_ELI_no_N_LFQ[,c(1:8)])),max(na.omit(YEAST_ELI_no_N_LFQ[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(YEAST_ELI_no_N_LFQ[,c(1:8)])),method="average")

pdf(file="QC_PME12_LFQ_YEAST_May19.pdf", width = 15, height = 15)
boxplot(YEAST_ELI_no_N_LFQ[,c(1:8)],names=colnames(YEAST_ELI_no_N_LFQ)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(YEAST_ELI_no_N_LFQ[,c(1:8)]+1), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_YEAST)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(YEAST_ELI_no_N_LFQ[,c(1:8)], tipo, "DetOut_PME12_LFQ_YEAST_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(YEAST_ELI_no_N_LFQ[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(YEAST_ELI_no_N_LFQ[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_LFQ_YEAST_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

################################ lo mismo para iBAQ

file<-read.csv2("/home/margaret/data/pepe/10_PME12_SPECTRAL_COUNTING_8_5_19/Results_MaxQuant/iBAQ_24_5_2019_MaxQuant.txt", header=T, sep="\t")


File_arr<-file[-c(1,2),]


#################### processing Perseus output

protein_groups<-unlist(File_arr$Protein.IDs)

only_1_specie_in_each_protein_group<-c()

for (i in 1:length(protein_groups)){
    print(i)
    char_b<-as.character(protein_groups[i])
    char_splitted<-unlist(strsplit(char_b,";"))
    protein<-paste(char_splitted[1])
    specie_to_search<-sapply(strsplit(protein,"_"),"[", 2)
    start_specie<-substr(specie_to_search, 1, 3)
    species_to_merge<-c()
    for( j in 1:length(char_splitted)){
    specie_to_merge<-sapply(strsplit(char_splitted[j],"_"),"[", 2)
    species_to_merge<-c(species_to_merge,specie_to_merge)
    }
    specie_to_merge_start<-c()
    for(k in 1:length(species_to_merge)){
        specie_to_merge_start[k]<-substr(species_to_merge[k], 1, 3)
        if(start_specie==specie_to_merge_start[k]){

            only_1_specie_in_each_protein_group[i]<-"TRUE"

        }else{
            only_1_specie_in_each_protein_group[i]<-"FALSE"

        }
    }
}

#### quitar groous con más de 1 specie

File_arr$only_one_specie<-paste(only_1_specie_in_each_protein_group)


data_final_selected<-File_arr[which(File_arr$only_one_specie=="TRUE"),]

specie_final<-c()
for (i in 1:nrow(data_final_selected)){
    print(i)
    char_b<-paste(as.character(data_final_selected$Protein.IDs[i]))
    char_splitted<-strsplit(char_b,";")
    protein<-char_splitted[[1]][1]
    specie_detected<-sapply(strsplit(protein,"_"),"[", 2)
    specie_final<-c(specie_final,specie_detected)
}
data_final_selected$ESPECIE<-paste(specie_final)




data_final_selected[,1]<-as.numeric(paste(data_final_selected[,1]))
data_final_selected[,2]<-as.numeric(paste(data_final_selected[,2]))
data_final_selected[,3]<-as.numeric(paste(data_final_selected[,3]))
data_final_selected[,4]<-as.numeric(paste(data_final_selected[,4]))
data_final_selected[,5]<-as.numeric(paste(data_final_selected[,5]))
data_final_selected[,6]<-as.numeric(paste(data_final_selected[,6]))
data_final_selected[,7]<-as.numeric(paste(data_final_selected[,7]))
data_final_selected[,8]<-as.numeric(paste(data_final_selected[,8]))





data_final_selected$media_iBAQ_A<- rowMeans(data_final_selected[,1:4])
data_final_selected$media_iBAQ_B<-rowMeans(data_final_selected[,5:8])

data_final_selected$RATIO_log2<-(data_final_selected$media_iBAQ_B-data_final_selected$media_iBAQ_A)

save(data_final_selected,file="quantification_Ratios_MAXQUANT_iBAQ_at_PSM_level.Rdata")


################ disgregar el objecto por especies

HUMAN<-data_final_selected[which(data_final_selected$ESPECIE =="HUMAN"),]   #3633
ECO<-data_final_selected[ grepl("^ECO",data_final_selected$ESPECIE),] #259
YEAST<-data_final_selected[grepl("^YE",data_final_selected$ESPECIE),] #1343
YEASTV<-data_final_selected[grepl("^S",data_final_selected$ESPECIE),]

YEAST<-rbind(YEAST,YEASTV)


proteins_good_quantified_HUMAN<-HUMAN[which((HUMAN$RATIO_log2 >=-0.514) & (HUMAN$RATIO_log2 <=0.378)),] #2347

proteins_good_quantified_YEAST<-YEAST[which((YEAST$RATIO_log2 >=0.485) & (YEAST$RATIO_log2 <=1.378)),] #632

proteins_good_quantified_ECO<-ECO[which((ECO$RATIO_log2 >= -2.514) &(ECO$RATIO_log2 <=-1.621 )),] #59

##### desviacion estandar  y dispersion

HUMAN$RATIO_RAW<-2^HUMAN$RATIO_log2
ECO$RATIO_RAW<-2^ECO$RATIO_log2
YEAST$RATIO_RAW<-2^YEAST$RATIO_log2

sd_HUMAN<-sd(HUMAN$RATIO_RAW)
mean_HUMAN<-mean(HUMAN$RATIO_RAW)

HUMAN_disp<-(sd_HUMAN/mean_HUMAN)*100

sd_ECO<-sd(ECO$RATIO_RAW)
mean_ECO<-mean(ECO$RATIO_RAW)

ECO_disp<-(sd_ECO/mean_ECO)*100

sd_YEAST<-sd(YEAST$RATIO_RAW)
mean_YEAST<-mean(YEAST$RATIO_RAW)

YEAST_disp<-(sd_YEAST/mean_YEAST)*100

#PLOTS


pdf("Distr_ratios_by_species_iBAQ.pdf")
par(mfrow=c(2,2))
hist(HUMAN$RATIO_RAW)
abline(v=c(0.7,1,1.3), col=c("red","blue","red"), lwd=2, lty=1)
hist(YEAST$RATIO_RAW)
abline(v=c(1.4,2,2.6), col=c("red","blue","red"), lwd=2, lty=1)
hist(ECO$RATIO_RAW)
abline(v=c(0.175,0.25,0.325), col=c("red","blue","red"), lwd=2, lty=1)
dev.off()

pdf("Distr_ratios_log2_by_species_iBAQ.pdf")
par(mfrow=c(2,2))
hist(HUMAN$RATIO_log2)
abline(v=c(-0.514,0,0.378), col=c("red","blue","red"), lwd=2, lty=1)
hist(YEAST$RATIO_log2)
abline(v=c(0.485,1,1.378), col=c("red","blue","red"), lwd=2, lty=1)
hist(ECO$RATIO_log2)
abline(v=c(-2.514,-2,-1.621), col=c("red","blue","red"), lwd=2, lty=1)
dev.off()

tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio_log"=HUMAN$RATIO_log2)
tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio_log"=ECO$RATIO_log2)
tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST)),"ratio_log"=YEAST$RATIO_log2)



plotter<-data.frame("organism"=NULL,"ratio_log"=NULL)
plotter<-rbind(plotter,tmp_HUMAN)
plotter<-rbind(plotter,tmp_ECO)
plotter<-rbind(plotter,tmp_YEAST)

#box `plot`
pdf("box_plot_distr_ratios_log_by_species_iBAQ.pdf")
ggplot(data=plotter, aes(x=organism, y=ratio_log, fill=organism)) + geom_boxplot() + geom_hline(yintercept=c(-2.514,-2,-1.621,-0.514,0,0.378,0.485,1,1.378), color=c("green","green","green","red","red","red","blue","blue","blue"))
dev.off()

# tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio"=HUMAN$RATIO_RAW)
# tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio"=ECO$RATIO_RAW)
# tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST)),"ratio"=YEAST$RATIO_RAW)
#
# library(ggplot2)
#
# plotter<-data.frame("organism"=NULL,"ratio"=NULL)
# plotter<-rbind(plotter,tmp_HUMAN)
# plotter<-rbind(plotter,tmp_ECO)
# plotter<-rbind(plotter,tmp_YEAST)
#
# #box `plot`
# pdf("box_plot_distr_ratios_by_species_iBAQ.pdf")
# ggplot(data=plotter, aes(x=organism, y=ratio, fill=organism)) + geom_boxplot() + geom_hline(yintercept=c(0.25,1,2), color=c("green","red","blue"))
# dev.off()

#scatter_plot

tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio_log"=HUMAN$RATIO_log2, "media_iBAQ_B"=HUMAN$media_iBAQ_B)
tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio_log"=ECO$RATIO_log2, "media_iBAQ_B"=ECO$media_iBAQ_B)
tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST)),"ratio_log"=YEAST$RATIO_log2, "media_iBAQ_B"=YEAST$media_iBAQ_B)



plotter<-data.frame("organism"=NULL,"ratio_log"=NULL, "media_iBAQ_B"=NULL)
plotter<-rbind(plotter,tmp_HUMAN)
plotter<-rbind(plotter,tmp_ECO)
plotter<-rbind(plotter,tmp_YEAST)



#scatter `plot`
pdf("scatter_plot_distr_ratios_log2_by_species_iBAQ.pdf")
ggplot(data=plotter, aes(x=media_iBAQ_B, y=ratio_log, fill=organism, color=organism)) + geom_point() + geom_hline(yintercept=c(-2.514,-2,-1.621,-0.514,0,0.378,0.485,1,1.378), color=c("green","green","green","red","red","red","blue","blue","blue"))
dev.off()

################333

#PLOTS DETOUTLIERS

#from data_final_selected para seleccionar quantificaciones basadas en 1 peptido

#all_one_peptide<-data_final_selected[which(data_final_selected$Peptides==1),]

# HUMAN<-all_one_peptide[which(all_one_peptide$ESPECIE =="HUMAN"),]   #3633
# ECO<-all_one_peptide[ grepl("^ECO",all_one_peptide$ESPECIE),] #259
# YEAST<-all_one_peptide[grepl("^YE",all_one_peptide$ESPECIE),] #1343
# YEASTV<-all_one_peptide[grepl("^S",all_one_peptide$ESPECIE),]
#
# YEAST<-rbind(YEAST,YEASTV)

#two_or_more_peptides<-data_final_selected[which(as.numeric(paste(data_final_selected$Peptides))>=2),]

HUMAN<-two_or_more_peptides[which(two_or_more_peptides$ESPECIE =="HUMAN"),]   #3633
ECO<-two_or_more_peptides[ grepl("^ECO",two_or_more_peptides$ESPECIE),] #259
YEAST<-two_or_more_peptides[grepl("^YE",two_or_more_peptides$ESPECIE),] #1343
YEASTV<-two_or_more_peptides[grepl("^S",two_or_more_peptides$ESPECIE),]

YEAST<-rbind(YEAST,YEASTV)


subset_ELI_no_N_iBAQ<-data_final_selected[,c(1:8,23)]
#subset_ELI_no_N_iBAQ<-all_one_peptide[,c(1:8,23)]
#subset_ELI_no_N_iBAQ<-two_or_more_peptides[,c(1:8,23)]





acol = sample(brewer.pal(8, "Dark2"), ncol(subset_ELI_no_N_iBAQ)-1, replace = (8 < ncol(subset_ELI_no_N_iBAQ)-1))
xlim = c(min(na.omit(subset_ELI_no_N_iBAQ[,c(1:8)])),max(na.omit(subset_ELI_no_N_iBAQ[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(subset_ELI_no_N_iBAQ[,c(1:8)])),method="average")

pdf(file="QC_PME12_iBAQ_May19.pdf", width = 15, height = 15)
boxplot(subset_ELI_no_N_iBAQ[,c(1:8)],names=colnames(subset_ELI_no_N_iBAQ)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(subset_ELI_no_N_iBAQ[,c(1:8)]+1), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_no_N_iBAQ[,c(1:8)], tipo, "DetOut_PME12_iBAQ_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_no_N_iBAQ[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_no_N_iBAQ[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_iBAQ_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

#### epecies separadas

HUMAN_ELI_no_N_iBAQ<-HUMAN[,c(1:8,23)]
ECO_ELI_no_N_iBAQ<-ECO[,c(1:8,23)]
YEAST_ELI_no_N_iBAQ<-YEAST[,c(1:8,23)]

####### HUMAN

acol = sample(brewer.pal(8, "Dark2"), ncol(HUMAN_ELI_no_N_iBAQ)-1, replace = (8 < ncol(HUMAN_ELI_no_N_iBAQ)-1))
xlim = c(min(na.omit(HUMAN_ELI_no_N_iBAQ[,c(1:8)])),max(na.omit(HUMAN_ELI_no_N_iBAQ[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(HUMAN_ELI_no_N_iBAQ[,c(1:8)])),method="average")

pdf(file="QC_PME12_iBAQ_HUMAN_May19.pdf", width = 15, height = 15)
boxplot(HUMAN_ELI_no_N_iBAQ[,c(1:8)],names=colnames(HUMAN_ELI_no_N_iBAQ)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(HUMAN_ELI_no_N_iBAQ[,c(1:8)]+1), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(HUMAN_ELI_no_N_iBAQ[,c(1:8)], tipo, "DetOut_PME12_iBAQ_HUMAN_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(HUMAN_ELI_no_N_iBAQ[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(HUMAN_ELI_no_N_iBAQ[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_iBAQ_HUMAN_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

#### ECOLI

acol = sample(brewer.pal(8, "Dark2"), ncol(ECO_ELI_no_N_iBAQ)-1, replace = (8 < ncol(ECO_ELI_no_N_iBAQ)-1))
xlim = c(min(na.omit(ECO_ELI_no_N_iBAQ[,c(1:8)])),max(na.omit(ECO_ELI_no_N_iBAQ[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(ECO_ELI_no_N_iBAQ[,c(1:8)])),method="average")

pdf(file="QC_PME12_iBAQ_ECO_May19.pdf", width = 15, height = 15)
boxplot(ECO_ELI_no_N_iBAQ[,c(1:8)],names=colnames(ECO_ELI_no_N_iBAQ)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(ECO_ELI_no_N_iBAQ[,c(1:8)]+1), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(ECO_ELI_no_N_iBAQ[,c(1:8)], tipo, "DetOut_PME12_iBAQ_ECO_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(ECO_ELI_no_N_iBAQ[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(ECO_ELI_no_N_iBAQ[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_iBAQ_ECO_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

### YEAST_

acol = sample(brewer.pal(8, "Dark2"), ncol(YEAST_ELI_no_N_iBAQ)-1, replace = (8 < ncol(YEAST_ELI_no_N_iBAQ)-1))
xlim = c(min(na.omit(YEAST_ELI_no_N_iBAQ[,c(1:8)])),max(na.omit(YEAST_ELI_no_N_iBAQ[,c(1:8)])))
clust.euclid.average <- hclust(dist(t(YEAST_ELI_no_N_iBAQ[,c(1:8)])),method="average")

pdf(file="QC_PME12_iBAQ_YEAST_May19.pdf", width = 15, height = 15)
boxplot(YEAST_ELI_no_N_iBAQ[,c(1:8)],names=colnames(YEAST_ELI_no_N_iBAQ)[c(1:8)], which="all", cex.axis=0.7, las=2)
multi("density",(YEAST_ELI_no_N_iBAQ[,c(1:8)]+1), xlim, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_YEAST)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(YEAST_ELI_no_N_iBAQ[,c(1:8)], tipo, "DetOut_PME12_iBAQ_YEAST_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(YEAST_ELI_no_N_iBAQ[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(YEAST_ELI_no_N_iBAQ[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_iBAQ_YEAST_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()
