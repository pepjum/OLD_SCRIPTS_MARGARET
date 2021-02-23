#WEIGHTED

results<-read.csv2("Protein_sets_union_21_05.csv", header=T, sep=",")
#"XIC_QUANTI_PEPTIDES.csv"

#results[,9:36]<-NULL
results_raw<-results[,c(2,37:60,29:36)]
results_raw$AC<-sapply(strsplit(paste(results_raw$description),"\\|"),"[", 2)
results_raw$ProteinName<-sapply(strsplit(paste(results_raw$description),"\\|"),"[", 3)
results_raw$ESPECIE<-sapply(strsplit(paste(results_raw$description),"_"),"[", 2)

results_raw_basic<-results_raw[,c("accession","Basic.SC_F304990","Basic.SC_F304991","Basic.SC_F304992","Basic.SC_F304993","Basic.SC_F304994","Basic.SC_F304995","Basic.SC_F304996","Basic.SC_F304997")]
write.table(results_raw_basic, "results_raw_basic.tsv", col.names=T, row.names=F, quote=F, sep="\t")

results_raw_Specific<-results_raw[,c("accession","Specific.SC_F304990","Specific.SC_F304991","Specific.SC_F304992","Specific.SC_F304993","Specific.SC_F304994","Specific.SC_F304995","Specific.SC_F304996","Specific.SC_F304997")]
write.table(results_raw_Specific, "results_raw_Specific.tsv", col.names=T, row.names=F, quote=F, sep="\t")

results_raw_Weighted<-results_raw[,c("accession","Weighted.SC_F304990","Weighted.SC_F304991","Weighted.SC_F304992","Weighted.SC_F304993","Weighted.SC_F304994","Weighted.SC_F304995","Weighted.SC_F304996","Weighted.SC_F304997")]
write.table(results_raw_Weighted, "results_raw_Weighted.tsv", col.names=T, row.names=F, quote=F, sep="\t")

# data_BASIC_SC<-results_raw[,c(34:36,2:9,26:33)]
# data_SPECIFIC_SC<-results_raw[,c(34:36,10:17,26:33)]
# data_WEIGHTED_SC<-results_raw[,c(34:36,18:25,26:33)]


subset_data_basic<-results_raw_basic
subset_data_specific<-results_raw_Specific
subset_data_weighted<-results_raw_Weighted

#eliminar las proteinas que no están en al menos 1 pinchazo en A y 1 en B


#results_sel<-results[,c("Protein.Set","Abundance.F304990","Abundance.F304991","Abundance.F304992","Abundance.F304993","Abundance.F304994","Abundance.F304995","Abundance.F304996","Abundance.F304997","Pep..match.count.F304990","Pep..match.count.F304991","Pep..match.count.F304992","Pep..match.count.F304993","Pep..match.count.F304994","Pep..match.count.F304995","Pep..match.count.F304996","Pep..match.count.F304997")]

##### BASIC SC eliminar las que no esten en al menos 1 cuadruplicado
#names(results_sel)[2:9]<-c("A_1","A_2","A_3","A_4","B_1","B_2","B_3","B_4")
# names(results_sel)[1]<-"ProteinAccession"
# results_sel$ESPECIE<-sapply(strsplit(paste(results_sel$ProteinAccession),"_"),"[", 2)
# results_sel$AC<-sapply(strsplit(paste(results_sel$ProteinAccession),"\\|"),"[", 2)
# results_sel$ProteinName<-sapply(strsplit(paste(results_sel$ProteinAccession),"\\|"),"[", 3)
#
# results_raw<-results_sel[,c("AC","ProteinName","ESPECIE","A_1","A_2","A_3","A_4","B_1","B_2","B_3","B_4")]
# write.table(results_raw, "results_raw.tsv", col.names=T,row.names=F, quote=F, sep="\t")


# results_sel<-results_sel[which(sum(results_sel[,2:9])>0),]
# results_sel$suma<-rowSums(results_sel[10:17])
# results_sel<-results_sel[which(results_sel$suma>=1),]
matriz_BASIC<-matrix(NA, length(unique(subset_data_basic$accession)), 8)
colnames(matriz_BASIC) <-names(subset_data_basic)[2:9]
rownames(matriz_BASIC) <- paste(unique(subset_data_basic$accession))

for(protein in unique(subset_data_basic$accession)){
    print(protein)
    tmp <- subset_data_basic[subset_data_basic$accession == protein, ]
    tmp<-unique(tmp)
    for (i in names(tmp[,2:9])){
    matriz_BASIC[protein, i] <- tmp[,i ]

    }
}
#
matriz_B<-matrix(NA, length(unique(subset_data_basic$accession)), 2)
matriz_all_BASIC<-cbind(matriz_BASIC,matriz_B)
colnames(matriz_all_BASIC)<-c(colnames(matriz_BASIC),"n_of_samples_quantified_A","n_of_samples_quantified_B")

for(protein in unique(subset_data_basic$accession)){
    print(protein)
    matriz_all_BASIC[protein,"n_of_samples_quantified_A"]<-sum(!is.na(matriz_all_BASIC[protein,][1:4]))
    matriz_all_BASIC[protein,"n_of_samples_quantified_B"]<-sum(!is.na(matriz_all_BASIC[protein,][5:8]))
}

save(matriz_all_BASIC, file="matriz_all_BASIC_SC_union_without_filtering_without_imputation.rda")

matriz_specific<-matrix(NA, length(unique(subset_data_specific$accession)), 8)
colnames(matriz_specific) <-names(subset_data_specific)[2:9]
rownames(matriz_specific) <- paste(unique(subset_data_specific$accession))

for(protein in unique(subset_data_specific$accession)){
    print(protein)
    tmp <- subset_data_specific[subset_data_specific$accession == protein, ]
    tmp<-unique(tmp)
    for (i in names(tmp[,2:9])){
    matriz_specific[protein, i] <- tmp[,i ]

    }
}
#
matriz_B<-matrix(NA, length(unique(subset_data_specific$accession)), 2)
matriz_all_specific<-cbind(matriz_specific,matriz_B)
colnames(matriz_all_specific)<-c(colnames(matriz_specific),"n_of_samples_quantified_A","n_of_samples_quantified_B")

for(protein in unique(subset_data_specific$accession)){
    print(protein)
    matriz_all_specific[protein,"n_of_samples_quantified_A"]<-sum(!is.na(matriz_all_specific[protein,][1:4]))
    matriz_all_specific[protein,"n_of_samples_quantified_B"]<-sum(!is.na(matriz_all_specific[protein,][5:8]))
}

save(matriz_all_specific, file="matriz_all_specific_SC_union_without_filtering_without_imputation.rda")


matriz_weighted<-matrix(NA, length(unique(subset_data_weighted$accession)), 8)
colnames(matriz_weighted) <-names(subset_data_weighted)[2:9]
rownames(matriz_weighted) <- paste(unique(subset_data_weighted$accession))

for(protein in unique(subset_data_weighted$accession)){
    print(protein)
    tmp <- subset_data_weighted[subset_data_weighted$accession == protein, ]
    tmp<-unique(tmp)
    for (i in names(tmp[,2:9])){
    matriz_weighted[protein, i] <- tmp[,i ]

    }
}
#
matriz_B<-matrix(NA, length(unique(subset_data_weighted$accession)), 2)
matriz_all_weighted<-cbind(matriz_weighted,matriz_B)
colnames(matriz_all_weighted)<-c(colnames(matriz_weighted),"n_of_samples_quantified_A","n_of_samples_quantified_B")

for(protein in unique(subset_data_weighted$accession)){
    print(protein)
    matriz_all_weighted[protein,"n_of_samples_quantified_A"]<-sum(!is.na(matriz_all_weighted[protein,][1:4]))
    matriz_all_weighted[protein,"n_of_samples_quantified_B"]<-sum(!is.na(matriz_all_weighted[protein,][5:8]))
}

save(matriz_all_weighted, file="matriz_all_weighted_SC_union_without_filtering_without_imputation.rda")



########## filtrar LOS QUE NO SE VEN AL MENOS 1 VEZ EN A Y 1 EN B

df_all_BASIC<-as.data.frame(matriz_all_BASIC)
df_all_BASIC$ProteinName<-rownames(matriz_all_BASIC)

df_all_filtered_by_nsamples_quantified_BASIC<-df_all_BASIC[which((df_all_BASIC$n_of_samples_quantified_A >=1) & (df_all_BASIC$n_of_samples_quantified_B >=1)  ),]
#
dat_BASIC<-df_all_filtered_by_nsamples_quantified_BASIC



df_all_SPECIFIC<-as.data.frame(matriz_all_specific)
df_all_SPECIFIC$ProteinName<-rownames(matriz_all_specific)

df_all_filtered_by_nsamples_quantified_SPECIFIC<-df_all_SPECIFIC[which((df_all_SPECIFIC$n_of_samples_quantified_A >=1) & (df_all_SPECIFIC$n_of_samples_quantified_B >=1)  ),]

dat_SPECIFIC<-df_all_filtered_by_nsamples_quantified_SPECIFIC


df_all_WEIGHTED<-as.data.frame(matriz_all_weighted)
df_all_WEIGHTED$ProteinName<-rownames(matriz_all_weighted)

df_all_filtered_by_nsamples_quantified_WEIGHTED<-df_all_WEIGHTED[which((df_all_WEIGHTED$n_of_samples_quantified_A >=1) & (df_all_WEIGHTED$n_of_samples_quantified_B >=1)  ),]

dat_WEIGHTED<-df_all_filtered_by_nsamples_quantified_WEIGHTED


save(dat_BASIC, file="matriz_all_BASIC_SC_union_filtered_without_imputation.rda")
save(dat_SPECIFIC, file="matriz_all_SPECIFIC_SC_union_filtered_without_imputation.rda")
save(dat_WEIGHTED, file="matriz_all_WEIGHTED_SC_union_filtered_without_imputation.rda")

# ###### NORMALIZAR POR PINCHAZO

dat_BASIC$normalized_A_1<-dat_BASIC[,1]-colMeans(dat_BASIC[1], na.rm=TRUE)
dat_BASIC$normalized_A_2<-dat_BASIC[,2]-colMeans(dat_BASIC[2], na.rm=TRUE)
dat_BASIC$normalized_A_3<-dat_BASIC[,3]-colMeans(dat_BASIC[3], na.rm=TRUE)
dat_BASIC$normalized_A_4<-dat_BASIC[,4]-colMeans(dat_BASIC[4], na.rm=TRUE)
dat_BASIC$normalized_B_1<-dat_BASIC[,5]-colMeans(dat_BASIC[5], na.rm=TRUE)
dat_BASIC$normalized_B_2<-dat_BASIC[,6]-colMeans(dat_BASIC[6], na.rm=TRUE)
dat_BASIC$normalized_B_3<-dat_BASIC[,7]-colMeans(dat_BASIC[7], na.rm=TRUE)
dat_BASIC$normalized_B_4<-dat_BASIC[,8]-colMeans(dat_BASIC[8], na.rm=TRUE)

dat_SPECIFIC$normalized_A_1<-dat_SPECIFIC[,1]-colMeans(dat_SPECIFIC[1], na.rm=TRUE)
dat_SPECIFIC$normalized_A_2<-dat_SPECIFIC[,2]-colMeans(dat_SPECIFIC[2], na.rm=TRUE)
dat_SPECIFIC$normalized_A_3<-dat_SPECIFIC[,3]-colMeans(dat_SPECIFIC[3], na.rm=TRUE)
dat_SPECIFIC$normalized_A_4<-dat_SPECIFIC[,4]-colMeans(dat_SPECIFIC[4], na.rm=TRUE)
dat_SPECIFIC$normalized_B_1<-dat_SPECIFIC[,5]-colMeans(dat_SPECIFIC[5], na.rm=TRUE)
dat_SPECIFIC$normalized_B_2<-dat_SPECIFIC[,6]-colMeans(dat_SPECIFIC[6], na.rm=TRUE)
dat_SPECIFIC$normalized_B_3<-dat_SPECIFIC[,7]-colMeans(dat_SPECIFIC[7], na.rm=TRUE)
dat_SPECIFIC$normalized_B_4<-dat_SPECIFIC[,8]-colMeans(dat_SPECIFIC[8], na.rm=TRUE)

dat_WEIGHTED$normalized_A_1<-dat_WEIGHTED[,1]-colMeans(dat_WEIGHTED[1], na.rm=TRUE)
dat_WEIGHTED$normalized_A_2<-dat_WEIGHTED[,2]-colMeans(dat_WEIGHTED[2], na.rm=TRUE)
dat_WEIGHTED$normalized_A_3<-dat_WEIGHTED[,3]-colMeans(dat_WEIGHTED[3], na.rm=TRUE)
dat_WEIGHTED$normalized_A_4<-dat_WEIGHTED[,4]-colMeans(dat_WEIGHTED[4], na.rm=TRUE)
dat_WEIGHTED$normalized_B_1<-dat_WEIGHTED[,5]-colMeans(dat_WEIGHTED[5], na.rm=TRUE)
dat_WEIGHTED$normalized_B_2<-dat_WEIGHTED[,6]-colMeans(dat_WEIGHTED[6], na.rm=TRUE)
dat_WEIGHTED$normalized_B_3<-dat_WEIGHTED[,7]-colMeans(dat_WEIGHTED[7], na.rm=TRUE)
dat_WEIGHTED$normalized_B_4<-dat_WEIGHTED[,8]-colMeans(dat_WEIGHTED[8], na.rm=TRUE)
#
# ######## MEDIA Y RATIO DE LOS VALORES NORMALIZADOS
#
dat_BASIC<-transform(dat_BASIC, media_BASIC_A_NORM = rowMeans(dat_BASIC[,12:15], na.rm = TRUE))
dat_BASIC<-transform(dat_BASIC, media_BASIC_B_NORM = rowMeans(dat_BASIC[,16:19], na.rm = TRUE))
dat_BASIC$RATIO_NORM<-dat_BASIC$media_BASIC_B_NORM/dat_BASIC$media_BASIC_A_NORM

dat_SPECIFIC<-transform(dat_SPECIFIC, media_SPECIFIC_A_NORM = rowMeans(dat_SPECIFIC[,12:15], na.rm = TRUE))
dat_SPECIFIC<-transform(dat_SPECIFIC, media_SPECIFIC_B_NORM = rowMeans(dat_SPECIFIC[,16:19], na.rm = TRUE))
dat_SPECIFIC$RATIO_NORM<-dat_SPECIFIC$media_SPECIFIC_B_NORM/dat_SPECIFIC$media_SPECIFIC_A_NORM

dat_WEIGHTED<-transform(dat_WEIGHTED, media_WEIGHTED_A_NORM = rowMeans(dat_WEIGHTED[,12:15], na.rm = TRUE))
dat_WEIGHTED<-transform(dat_WEIGHTED, media_WEIGHTED_B_NORM = rowMeans(dat_WEIGHTED[,16:19], na.rm = TRUE))
dat_WEIGHTED$RATIO_NORM<-dat_WEIGHTED$media_WEIGHTED_B_NORM/dat_WEIGHTED$media_WEIGHTED_A_NORM

save(dat_BASIC, file="matriz_all_BASIC_SC_filtered_without_imputation_RATIOS_NORMALIZED.rda")
save(dat_SPECIFIC, file="matriz_all_SPECIFIC_SC_filtered_without_imputation_RATIOS_NORMALIZED.rda")
save(dat_WEIGHTED, file="matriz_all_WEIGHTED_SC_filtered_without_imputation_RATIOS_NORMALIZED.rda")
#
#
# #dat_BASIC$RATIO_NO_NORM[which(dat_BASIC$RATIO_NO_NORM ==Inf)]<-10.8333333333333
#
# media_8_pinchazos_BASIC<-mean(colSums(dat_BASIC[1:8], na.rm=TRUE))
#
# dat_BASIC$normalized_A_1_ELI<-dat_BASIC[,1]/ (sum(dat_BASIC[,1],na.rm=TRUE)/media_8_pinchazos_BASIC)
# dat_BASIC$normalized_A_2_ELI<-dat_BASIC[,2]/ (sum(dat_BASIC[,2],na.rm=TRUE)/media_8_pinchazos_BASIC)
# dat_BASIC$normalized_A_3_ELI<-dat_BASIC[,3]/ (sum(dat_BASIC[,3],na.rm=TRUE)/media_8_pinchazos_BASIC)
# dat_BASIC$normalized_A_4_ELI<-dat_BASIC[,4]/ (sum(dat_BASIC[,4],na.rm=TRUE)/media_8_pinchazos_BASIC)
# dat_BASIC$normalized_B_1_ELI<-dat_BASIC[,5]/ (sum(dat_BASIC[,5],na.rm=TRUE)/media_8_pinchazos_BASIC)
# dat_BASIC$normalized_B_2_ELI<-dat_BASIC[,6]/ (sum(dat_BASIC[,6],na.rm=TRUE)/media_8_pinchazos_BASIC)
# dat_BASIC$normalized_B_3_ELI<-dat_BASIC[,7]/ (sum(dat_BASIC[,7],na.rm=TRUE)/media_8_pinchazos_BASIC)
# dat_BASIC$normalized_B_4_ELI<-dat_BASIC[,8]/ (sum(dat_BASIC[,8],na.rm=TRUE)/media_8_pinchazos_BASIC)
#
#
# dat_BASIC<-transform(dat_BASIC, media_BASIC_A_ELINORM = rowMeans(dat_BASIC[,25:29], na.rm = TRUE))
# dat_BASIC<-transform(dat_BASIC, media_BASIC_B_ELINORM = rowMeans(dat_BASIC[,30:33], na.rm = TRUE))
# dat_BASIC$RATIO_ELINORM<-dat_BASIC$media_BASIC_B_ELINORM/dat_BASIC$media_BASIC_A_ELINORM
#
#
# media_8_pinchazos_SPECIFIC<-mean(colSums(dat_SPECIFIC[1:8], na.rm=TRUE))
#
# dat_SPECIFIC$normalized_A_1_ELI<-dat_SPECIFIC[,1]/ (sum(dat_SPECIFIC[,1],na.rm=TRUE)/media_8_pinchazos_SPECIFIC)
# dat_SPECIFIC$normalized_A_2_ELI<-dat_SPECIFIC[,2]/ (sum(dat_SPECIFIC[,2],na.rm=TRUE)/media_8_pinchazos_SPECIFIC)
# dat_SPECIFIC$normalized_A_3_ELI<-dat_SPECIFIC[,3]/ (sum(dat_SPECIFIC[,3],na.rm=TRUE)/media_8_pinchazos_SPECIFIC)
# dat_SPECIFIC$normalized_A_4_ELI<-dat_SPECIFIC[,4]/ (sum(dat_SPECIFIC[,4],na.rm=TRUE)/media_8_pinchazos_SPECIFIC)
# dat_SPECIFIC$normalized_B_1_ELI<-dat_SPECIFIC[,5]/ (sum(dat_SPECIFIC[,5],na.rm=TRUE)/media_8_pinchazos_SPECIFIC)
# dat_SPECIFIC$normalized_B_2_ELI<-dat_SPECIFIC[,6]/ (sum(dat_SPECIFIC[,6],na.rm=TRUE)/media_8_pinchazos_SPECIFIC)
# dat_SPECIFIC$normalized_B_3_ELI<-dat_SPECIFIC[,7]/ (sum(dat_SPECIFIC[,7],na.rm=TRUE)/media_8_pinchazos_SPECIFIC)
# dat_SPECIFIC$normalized_B_4_ELI<-dat_SPECIFIC[,8]/ (sum(dat_SPECIFIC[,8],na.rm=TRUE)/media_8_pinchazos_SPECIFIC)
#
#
# dat_SPECIFIC<-transform(dat_SPECIFIC, media_SPECIFIC_A_ELINORM = rowMeans(dat_SPECIFIC[,25:29], na.rm = TRUE))
# dat_SPECIFIC<-transform(dat_SPECIFIC, media_SPECIFIC_B_ELINORM = rowMeans(dat_SPECIFIC[,30:33], na.rm = TRUE))
# dat_SPECIFIC$RATIO_ELINORM<-dat_SPECIFIC$media_SPECIFIC_B_ELINORM/dat_SPECIFIC$media_SPECIFIC_A_ELINORM
#
#
# media_8_pinchazos_WEIGHTED<-mean(colSums(dat_WEIGHTED[1:8], na.rm=TRUE))
#
# dat_WEIGHTED$normalized_A_1_ELI<-dat_WEIGHTED[,1]/ (sum(dat_WEIGHTED[,1],na.rm=TRUE)/media_8_pinchazos_WEIGHTED)
# dat_WEIGHTED$normalized_A_2_ELI<-dat_WEIGHTED[,2]/ (sum(dat_WEIGHTED[,2],na.rm=TRUE)/media_8_pinchazos_WEIGHTED)
# dat_WEIGHTED$normalized_A_3_ELI<-dat_WEIGHTED[,3]/ (sum(dat_WEIGHTED[,3],na.rm=TRUE)/media_8_pinchazos_WEIGHTED)
# dat_WEIGHTED$normalized_A_4_ELI<-dat_WEIGHTED[,4]/ (sum(dat_WEIGHTED[,4],na.rm=TRUE)/media_8_pinchazos_WEIGHTED)
# dat_WEIGHTED$normalized_B_1_ELI<-dat_WEIGHTED[,5]/ (sum(dat_WEIGHTED[,5],na.rm=TRUE)/media_8_pinchazos_WEIGHTED)
# dat_WEIGHTED$normalized_B_2_ELI<-dat_WEIGHTED[,6]/ (sum(dat_WEIGHTED[,6],na.rm=TRUE)/media_8_pinchazos_WEIGHTED)
# dat_WEIGHTED$normalized_B_3_ELI<-dat_WEIGHTED[,7]/ (sum(dat_WEIGHTED[,7],na.rm=TRUE)/media_8_pinchazos_WEIGHTED)
# dat_WEIGHTED$normalized_B_4_ELI<-dat_WEIGHTED[,8]/ (sum(dat_WEIGHTED[,8],na.rm=TRUE)/media_8_pinchazos_WEIGHTED)
#
#
# dat_WEIGHTED<-transform(dat_WEIGHTED, media_WEIGHTED_A_ELINORM = rowMeans(dat_WEIGHTED[,25:29], na.rm = TRUE))
# dat_WEIGHTED<-transform(dat_WEIGHTED, media_WEIGHTED_B_ELINORM = rowMeans(dat_WEIGHTED[,30:33], na.rm = TRUE))
# dat_WEIGHTED$RATIO_ELINORM<-dat_WEIGHTED$media_WEIGHTED_B_ELINORM/dat_WEIGHTED$media_WEIGHTED_A_ELINORM
#
#
# save(dat_BASIC, file="matriz_all_BASIC_SC_filtered_without_imputation_RATIOS_NORMALIZED_AND_NON_NORMALIZED_AND_ELINORMALIZED.rda")
# save(dat_SPECIFIC, file="matriz_all_SPECIFIC_SC_filtered_without_imputation_RATIOS_NORMALIZED_AND_NON_NORMALIZED_AND_ELINORMALIZED.rda")
# save(dat_WEIGHTED, file="matriz_all_WEIGHTED_SC_filtered_without_imputation_RATIOS_NORMALIZED_AND_NON_NORMALIZED_AND_ELINORMALIZED.rda")

#save(results_sel,file="XIC_results_raw.rdata")

#medias XIC

# results_sel<-transform(results_sel, media_A = rowMeans(results_sel[,2:5]))
# results_sel<-transform(results_sel, media_B = rowMeans(results_sel[,6:9]))
#
# results_sel$log2_media_A<-log2(results_sel$media_A)
# results_sel$log2_media_B<-log2(results_sel$media_B)
#
# results_sel$ratio<-results_sel$log2_media_B-results_sel$log2_media_A
# names(results_sel)[1]<-"ProteinAccession"
# #######################################################
# ### detecto especie
# results_sel$ESPECIE<-sapply(strsplit(paste(results_sel$ProteinAccession),"_"),"[", 2)
# results_sel$AC<-sapply(strsplit(paste(results_sel$ProteinAccession),"\\|"),"[", 2)
# results_sel$ProteinName<-sapply(strsplit(paste(results_sel$ProteinAccession),"\\|"),"[", 3)



dat_BASIC$ESPECIE<-sapply(strsplit(paste(dat_BASIC$ProteinName),"_"),"[",2)
dat_SPECIFIC$ESPECIE<-sapply(strsplit(paste(dat_SPECIFIC$ProteinName),"_"), "[", 2)
dat_WEIGHTED$ESPECIE<-sapply(strsplit(paste(dat_WEIGHTED$ProteinName),"_"), "[", 2)
# write.table(results_sel, "processed.tsv", col.names=T, row.names=F, quote=F, sep="\t")
#separo BASIC por especies

# HUMAN<-results_sel[which(results_sel$ESPECIE =="HUMAN"),]   #3134
# ECO<-results_sel[ grepl("^ECO",results_sel$ESPECIE),] #199
# YEASTC<-results_sel[grepl("^YE",results_sel$ESPECIE),]
# YEASTV<-results_sel[grepl("^S",results_sel$ESPECIE),]
#
# YEAST<-rbind(YEASTC,YEASTV)  #1152

#proteins<-read.csv2("PROTEINSXIC.csv", header=T, sep=",")


HUMAN_BASIC<-dat_BASIC[which(dat_BASIC$ESPECIE =="HUMAN"),]   #3134
ECO_BASIC<-dat_BASIC[ grepl("^ECO",dat_BASIC$ESPECIE),] #199
YEAST_BASIC<-dat_BASIC[grepl("^YE",dat_BASIC$ESPECIE),]
YEASTV_BASIC<-dat_BASIC[grepl("^S",dat_BASIC$ESPECIE),]

YEAST_BASIC<-rbind(YEAST_BASIC,YEASTV_BASIC)  #1152


# proteins_good_quantified_HUMAN<-HUMAN[which((HUMAN$ratio >= -0.514) & (HUMAN$ratio <=0.378)),] #2431
# proteins_good_quantified_ECO<-ECO[which((ECO$ratio >= -2.514) & (ECO$ratio <=-1.621)),] #2431
# proteins_good_quantified_YEAST<-YEAST[which((YEAST$ratio >= 0.485) & (YEAST$ratio <=1.378)),] #2431
#
# single_pept_quant_good_HUMAN<-proteins_good_quantified_HUMAN[which((proteins_good_quantified_HUMAN$Pep..match.count.F304990 <=1)& (proteins_good_quantified_HUMAN$Pep..match.count.F304991 <=1) & (proteins_good_quantified_HUMAN$Pep..match.count.F304992 <=1) &(proteins_good_quantified_HUMAN$Pep..match.count.F304993 <=1) & (proteins_good_quantified_HUMAN$Pep..match.count.F304994 <=1) & (proteins_good_quantified_HUMAN$Pep..match.count.F304995 <=1) & (proteins_good_quantified_HUMAN$Pep..match.count.F304996 <=1) & (proteins_good_quantified_HUMAN$Pep..match.count.F304997 <=1) ),]
#
#
# single_pept_quant_good_ECO<-proteins_good_quantified_ECO[which((proteins_good_quantified_ECO$Pep..match.count.F304990 <=1)& (proteins_good_quantified_ECO$Pep..match.count.F304991 <=1) & (proteins_good_quantified_ECO$Pep..match.count.F304992 <=1) &(proteins_good_quantified_ECO$Pep..match.count.F304993 <=1) & (proteins_good_quantified_ECO$Pep..match.count.F304994 <=1) & (proteins_good_quantified_ECO$Pep..match.count.F304995 <=1) & (proteins_good_quantified_ECO$Pep..match.count.F304996 <=1) & (proteins_good_quantified_ECO$Pep..match.count.F304997 <=1) ),]
#
# single_pept_quant_good_YEAST<-proteins_good_quantified_YEAST[which((proteins_good_quantified_YEAST$Pep..match.count.F304990 <=1)& (proteins_good_quantified_YEAST$Pep..match.count.F304991 <=1) & (proteins_good_quantified_YEAST$Pep..match.count.F304992 <=1) &(proteins_good_quantified_YEAST$Pep..match.count.F304993 <=1) & (proteins_good_quantified_YEAST$Pep..match.count.F304994 <=1) & (proteins_good_quantified_YEAST$Pep..match.count.F304995 <=1) & (proteins_good_quantified_YEAST$Pep..match.count.F304996 <=1) & (proteins_good_quantified_YEAST$Pep..match.count.F304997 <=1) ),]


# sd
# HUMAN$RATIO_RAW<-2^HUMAN$ratio
# ECO$RATIO_RAW<-2^ECO$ratio
# YEAST$RATIO_RAW<-2^YEAST$ratio
# HUMAN$RATIO_RAW[is.infinite(HUMAN$RATIO_RAW)]<-NA
# ECO$RATIO_RAW[is.infinite(ECO$RATIO_RAW)]<-NA
# YEAST$RATIO_RAW[is.infinite(YEAST$RATIO_RAW)]<-NA
#
#
#
# sd_HUMAN<-sd(HUMAN$RATIO_RAW, na.rm=TRUE)
# mean_HUMAN<-mean(HUMAN$RATIO_RAW, na.rm=TRUE)
#
# HUMAN_disp<-(sd_HUMAN/mean_HUMAN)*100
#
#
# sd_ECO<-sd(ECO$RATIO_RAW, na.rm=TRUE)
# mean_ECO<-mean(ECO$RATIO_RAW, na.rm=TRUE)
#
# ECO_disp<-(sd_ECO/mean_ECO)*100
#
# sd_YEAST<-sd(YEAST$RATIO_RAW, na.rm=TRUE)
# mean_YEAST<-mean(YEAST$RATIO_RAW, na.rm=TRUE)
#
# YEAST_disp<-(sd_YEAST/mean_YEAST)*100

#plot(SPC_proteins_HUMAN$RATIO_log10)




proteins_good_quantified_HUMAN_BASIC<-HUMAN_BASIC[which((HUMAN_BASIC$RATIO_NORM >= 0.7) & (HUMAN_BASIC$RATIO_NORM <=1.3)),] #2431

proteins_good_quantified_YEAST_BASIC<-YEAST_BASIC[which((YEAST_BASIC$RATIO_NORM >=1.40) & (YEAST_BASIC$RATIO_NORM <=2.60)),] #3

proteins_good_quantified_ECO_BASIC<-ECO_BASIC[which((ECO_BASIC$RATIO_NORM >= 0.175) &(ECO_BASIC$RATIO_NORM <=0.325 )),] #4

#separo SPECIFIC POR especies

# HUMAN_SPECIFIC<-dat_SPECIFIC[which(dat_SPECIFIC$ESPECIE =="HUMAN"),]   #3134
# ECO_SPECIFIC<-dat_SPECIFIC[ grepl("^ECO",dat_SPECIFIC$ESPECIE),] #199
# YEAST_SPECIFIC<-dat_SPECIFIC[grepl("^YE",dat_SPECIFIC$ESPECIE),]
# YEASTV_SPECIFIC<-dat_SPECIFIC[grepl("^S",dat_SPECIFIC$ESPECIE),]
#
# YEAST_SPECIFIC<-rbind(YEAST_SPECIFIC,YEASTV_SPECIFIC)  #1152
#
# proteins_good_quantified_HUMAN_SPECIFIC<-HUMAN_SPECIFIC[which((HUMAN_SPECIFIC$RATIO_NO_NORM >= 0.7) & (HUMAN_SPECIFIC$RATIO_NO_NORM <=1.3)),] #633
#
# proteins_good_quantified_YEAST_SPECIFIC<-YEAST_SPECIFIC[which((YEAST_SPECIFIC$RATIO_NO_NORM >=1.40) & (YEAST_SPECIFIC$RATIO_NO_NORM <=2.60)),] #608
#
# proteins_good_quantified_ECO_SPECIFIC<-ECO_SPECIFIC[which((ECO_SPECIFIC$RATIO_NO_NORM >= 0.175) &(ECO_SPECIFIC$RATIO_NO_NORM <=0.325 )),] #10
#
# #separo WEIGHTED por ESPECIES
#
# HUMAN_WEIGHTED<-dat_WEIGHTED[which(dat_WEIGHTED$ESPECIE =="HUMAN"),]   #3134
# ECO_WEIGHTED<-dat_WEIGHTED[ grepl("^ECO",dat_WEIGHTED$ESPECIE),] #199
# YEAST_WEIGHTED<-dat_WEIGHTED[grepl("^YE",dat_WEIGHTED$ESPECIE),]
# YEASTV_WEIGHTED<-dat_WEIGHTED[grepl("^S",dat_WEIGHTED$ESPECIE),]
#
# YEAST_WEIGHTED<-rbind(YEAST_WEIGHTED,YEASTV_WEIGHTED)  #1152
#
# proteins_good_quantified_HUMAN_WEIGHTED<-HUMAN_WEIGHTED[which((HUMAN_WEIGHTED$RATIO_NO_NORM >= 0.7) & (HUMAN_WEIGHTED$RATIO_NO_NORM <=1.3)),] #679
#
# proteins_good_quantified_YEAST_WEIGHTED<-YEAST_WEIGHTED[which((YEAST_WEIGHTED$RATIO_NO_NORM >=1.40) & (YEAST_WEIGHTED$RATIO_NO_NORM <=2.60)),] #602
#
# proteins_good_quantified_ECO_WEIGHTED<-ECO_WEIGHTED[which((ECO_WEIGHTED$RATIO_NO_NORM >= 0.175) &(ECO_WEIGHTED$RATIO_NO_NORM <=0.325 )),] #7
#
#
# HUMAN_BASIC$RATIO_log2_NO_NORM<-log2(HUMAN_BASIC$RATIO_NO_NORM)
#
# HUMAN_SPECIFIC$RATIO_log2_NO_NORM<-log2(HUMAN_SPECIFIC$RATIO_NO_NORM)
# HUMAN_WEIGHTED$RATIO_log2_NO_NORM<-log2(HUMAN_WEIGHTED$RATIO_NO_NORM)
#
##############################################################################################################

####PLOTS

#########

# tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN_BASIC)),"ratio"=(HUMAN_BASIC$RATIO_ELINORM), "media_BASIC_B"=HUMAN_BASIC$media_BASIC_B_ELINORM)
# tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO_BASIC)),"ratio"=(ECO_BASIC$RATIO_ELINORM), "media_BASIC_B"=ECO_BASIC$media_BASIC_B_ELINORM)
# tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST_BASIC)),"ratio"=(YEAST_BASIC$RATIO_ELINORM), "media_BASIC_B"=YEAST_BASIC$media_BASIC_B_ELINORM)

pdf("Distr_ratios_by_species_XIC.pdf")
par(mfrow=c(2,2))
hist(HUMAN$RATIO_RAW)
abline(v=c(0.7,1,1.3), col=c("red","blue","red"), lwd=2, lty=1)
hist(YEAST$RATIO_RAW)
abline(v=c(1.4,2,2.6), col=c("red","blue","red"), lwd=2, lty=1)
hist(ECO$RATIO_RAW)
abline(v=c(0.175,0.25,0.325), col=c("red","blue","red"), lwd=2, lty=1)
dev.off()

pdf("Distr_ratios_log2_by_species_XIC.pdf")
par(mfrow=c(2,2))
hist(HUMAN$ratio)
abline(v=c(-0.514,0,0.378), col=c("red","blue","red"), lwd=2, lty=1)
hist(YEAST$ratio)
abline(v=c(0.485,1,1.378), col=c("red","blue","red"), lwd=2, lty=1)
hist(ECO$ratio)
abline(v=c(-2.514,-2,-1.621), col=c("red","blue","red"), lwd=2, lty=1)
dev.off()



















tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN)),"ratio"=(HUMAN$ratio), "media_B"=HUMAN$media_B)
tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO)),"ratio"=(ECO$ratio), "media_B"=ECO$media_B)
tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST)),"ratio"=(YEAST$ratio), "media_B"=YEAST$media_B)


library(ggplot2)

plotter_BASIC<-data.frame("organism"=NULL,"ratio"=NULL, "media_B"=NULL)
plotter_BASIC<-rbind(plotter_BASIC,tmp_HUMAN)
plotter_BASIC<-rbind(plotter_BASIC,tmp_ECO)
plotter_BASIC<-rbind(plotter_BASIC,tmp_YEAST)
plotter_BASIC$media_B<-log2(plotter_BASIC$media_B)

pdf("scatter_plot_distr_ratios_by_species_XIC.pdf")
par(mfrow=c(1,1))
ggplot(data=plotter_BASIC, aes(x=media_B, y=ratio, fill=organism, color=organism)) + geom_point() + geom_hline(yintercept=c(-2.514,-2,-1.621,-0.514,0,0.378,0.485,1,1.378), color=c("green","green","green","red","red","red","blue","blue","blue"), linetype=c("dashed","solid","dashed","dashed","solid","dashed","dashed","solid","dashed"))
dev.off()




# tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN_SPECIFIC)),"ratio"=(HUMAN_SPECIFIC$RATIO_ELINORM), "media_SPECIFIC_B"=HUMAN_SPECIFIC$media_SPECIFIC_B_ELINORM)
# tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO_SPECIFIC)),"ratio"=(ECO_SPECIFIC$RATIO_ELINORM), "media_SPECIFIC_B"=ECO_SPECIFIC$media_SPECIFIC_B_ELINORM)
# tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST_SPECIFIC)),"ratio"=(YEAST_SPECIFIC$RATIO_ELINORM), "media_SPECIFIC_B"=YEAST_SPECIFIC$media_SPECIFIC_B_ELINORM)
#
# plotter_SPECIFIC<-data.frame("organism"=NULL,"ratio"=NULL, "media_SPECIFIC_B"=NULL)
# plotter_SPECIFIC<-rbind(plotter_SPECIFIC,tmp_HUMAN)
# plotter_SPECIFIC<-rbind(plotter_SPECIFIC,tmp_ECO)
# plotter_SPECIFIC<-rbind(plotter_SPECIFIC,tmp_YEAST)
#
# tmp_HUMAN<-data.frame("organism"=rep("HeLa",nrow(HUMAN_WEIGHTED)),"ratio"=(HUMAN_WEIGHTED$RATIO_ELINORM), "media_WEIGHTED_B"=HUMAN_WEIGHTED$media_WEIGHTED_B_ELINORM)
# tmp_ECO<-data.frame("organism"=rep("E.coli",nrow(ECO_WEIGHTED)),"ratio"=(ECO_WEIGHTED$RATIO_ELINORM), "media_WEIGHTED_B"=ECO_WEIGHTED$media_WEIGHTED_B_ELINORM)
# tmp_YEAST<-data.frame("organism"=rep("Yeast",nrow(YEAST_WEIGHTED)),"ratio"=(YEAST_WEIGHTED$RATIO_ELINORM), "media_WEIGHTED_B"=YEAST_WEIGHTED$media_WEIGHTED_B_ELINORM)
#
# plotter_WEIGHTED<-data.frame("organism"=NULL,"ratio"=NULL, "media_WEIGHTED_B"=NULL)
# plotter_WEIGHTED<-rbind(plotter_WEIGHTED,tmp_HUMAN)
# plotter_WEIGHTED<-rbind(plotter_WEIGHTED,tmp_ECO)
# plotter_WEIGHTED<-rbind(plotter_WEIGHTED,tmp_YEAST)
#
#
#
#
# #scatter `plot`
# pdf("scatter_plot_distr_ratios_by_species_UNION_ELINORM.pdf")
# par(mfrow=c(3,1))
# ggplot(data=plotter_BASIC, aes(x=media_BASIC_B, y=ratio, fill=organism, color=organism)) + geom_point() + geom_hline(yintercept=c(0.175,0.25,0.325,0.7,1,1.3,1.4,2,2.6), color=c("green","green","green","red","red","red","blue","blue","blue"))
# ggplot(data=plotter_SPECIFIC, aes(x=media_SPECIFIC_B, y=ratio, fill=organism, color=organism)) + geom_point() + geom_hline(yintercept=c(0.175,0.25,0.325,0.7,1,1.3,1.4,2,2.6), color=c("green","green","green","red","red","red","blue","blue","blue"))
# ggplot(data=plotter_WEIGHTED, aes(x=media_WEIGHTED_B, y=ratio, fill=organism, color=organism)) + geom_point() + geom_hline(yintercept=c(0.175,0.25,0.325,0.7,1,1.3,1.4,2,2.6), color=c("green","green","green","red","red","red","blue","blue","blue"))
#
#
# dev.off()

# tmp_pairs_Basic<-dat_BASIC[,c(1:8)]
#
# library(psych)
# pdf("pairs_all.pdf")
# pairs.panels(tmp_pairs_Basic,
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = TRUE # show correlation ellipses
#              )
#
#
# dev.off()
#
# tmp_HUMAN_BASIC<-HUMAN_BASIC[,c(1:8)]
# tmp_ECO_BASIC<-ECO_BASIC[,c(1:8)]
# tmp_YEAST_BASIC<-YEAST_BASIC[,c(1:8)]
#
# pdf("pairs_by_species_BASIC.pdf")
# par(mfrow=c(1,3))
# pairs.panels(tmp_HUMAN_BASIC,
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = TRUE # show correlation ellipses
#              )
#
# pairs.panels(tmp_ECO_BASIC,
#             method = "pearson", # correlation method
#             hist.col = "#00AFBB",
#             density = TRUE,  # show density plots
#             ellipses = TRUE # show correlation ellipses
#             )
#
# pairs.panels(tmp_YEAST_BASIC,
#             method = "pearson", # correlation method
#             hist.col = "#00AFBB",
#             density = TRUE,  # show density plots
#             ellipses = TRUE # show correlation ellipses
#             )
#
#
# dev.off()
####### PLOTS DE ELI

#NORMALIZACION VS NO NORMALIZAR

#subset_ELI_no_N_BASIC<-dat_BASIC[,c(1:8,11)]
#subset_ELI_NORM_BASIC<-dat_BASIC[,c(26:33,11)]

subset_ELI_no_N_BASIC<-results_sel[,c(2:9,1)]

#subset_ELI_NORM_BASIC[is.na(subset_ELI_NORM_BASIC)]<-0


acol = sample(brewer.pal(8, "Dark2"), ncol(subset_ELI_no_N_BASIC)-1, replace = (8 < ncol(subset_ELI_no_N_BASIC)-1))
xlim_r = c(min(na.omit(log2(subset_ELI_no_N_BASIC[,c(1:8)]+1))),max(na.omit(log2(subset_ELI_no_N_BASIC[,c(1:8)]+1))))
xlim = c(min(na.omit(subset_ELI_no_N_BASIC[,c(1:8)]+1)),max(na.omit(subset_ELI_no_N_BASIC[,c(1:8)]+1)))
clust.euclid.average <- hclust(dist(t(log2(subset_ELI_no_N_BASIC[,c(1:8)]+1))),method="average")
#clust.euclid.average2 <- hclust(dist(t(na.omit(subset_ELI_NORM_BASIC[,c(1:8)]))),method="average")


#subset_ELI_no_N_BASIC<-dataMatCounts_f

source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesVikv2.R")

pdf(file="QC_PME12_May19.pdf", width = 15, height = 15)
boxplot(log2(subset_ELI_no_N_BASIC[,c(1:8)]+1),names=colnames(subset_ELI_no_N_BASIC)[c(1:8)], which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(subset_ELI_no_N_BASIC[,c(1:8)]+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
#boxplot(subset_ELI_NORM_BASIC[,c(1:8)], cex.axis=0.7, las=2)
#multi("density", subset_ELI_NORM_BASIC[,c(1:8)], xlim, "Density","","")
#plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_no_N_BASIC[,c(1:8)], tipo, "DetOut_PME12_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_no_N_BASIC[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_no_N_BASIC[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

#### lo mismo por especies
setwd("/home/margaret/data/pepe/10_PME12_SPECTRAL_COUNTING_8_5_19/QUALITY_CONTROLS/")

subset_ELI_no_N_BASIC_HUMAN<-HUMAN_BASIC[,c(1:8,11)]
subset_ELI_NORM_BASIC_HUMAN<-HUMAN_BASIC[,c(26:33,11)]

subset_ELI_no_N_BASIC_YEAST<-YEAST_BASIC[,c(1:8,11)]
subset_ELI_NORM_BASIC_YEAST<-YEAST_BASIC[,c(26:33,11)]

subset_ELI_no_N_BASIC_ECO<-ECO_BASIC[,c(1:8,11)]
subset_ELI_NORM_BASIC_ECO<-ECO_BASIC[,c(26:33,11)]

subset_ELI_no_N_HUMAN<-HUMAN[,c(2:9,1)]
subset_ELI_no_N_ECO<-ECO[,c(2:9,1)]
subset_ELI_no_N_YEAST<-YEAST[,c(2:9,1)]



#############
acol = sample(brewer.pal(8, "Dark2"), ncol(subset_ELI_no_N_HUMAN)-1, replace = (8 < ncol(subset_ELI_no_N_HUMAN)-1))
xlim_r = c(min(na.omit(log2(subset_ELI_no_N_HUMAN[,c(1:8)]+1))),max(na.omit(log2(subset_ELI_no_N_HUMAN[,c(1:8)]+1))))
xlim = c(min(na.omit(subset_ELI_no_N_HUMAN[,c(1:8)]+1)),max(na.omit(subset_ELI_no_N_HUMAN[,c(1:8)]+1)))

clust.euclid.average <- hclust(dist(t(log2(subset_ELI_no_N_HUMAN[,c(1:8)]+1))),method="average")
#clust.euclid.average2 <- hclust(dist(t(na.omit(subset_ELI_NORM_BASIC_HUMAN[,c(1:8)]))),method="average")


pdf(file="QC_PME12_HUMAN_May19.pdf", width = 15, height = 15)
boxplot(log2(subset_ELI_no_N_HUMAN[,c(1:8)]+1),names=colnames(subset_ELI_no_N_HUMAN)[c(1:8)], which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(subset_ELI_no_N_HUMAN[,c(1:8)]+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
#boxplot(log2(subset_ELI_NORM_BASIC_HUMAN[,c(1:8)]+1), cex.axis=0.7, las=2)
#multi("density", subset_ELI_NORM_BASIC_HUMAN[,c(1:8)], xlim, "Density","","")
#plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_HUMAN)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_no_N_HUMAN[,c(1:8)], tipo, "DetOut_PME12_HUMAN_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_no_N_HUMAN[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_no_N_HUMAN[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_HUMAN_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

acol = sample(brewer.pal(8, "Dark2"), ncol(subset_ELI_no_N_ECO)-1, replace = (8 < ncol(subset_ELI_no_N_ECO)-1))
xlim_r = c(min(na.omit(log2(subset_ELI_no_N_ECO[,c(1:8)]+1))),max(na.omit(log2(subset_ELI_no_N_ECO[,c(1:8)]+1))))
xlim = c(min(na.omit(subset_ELI_no_N_ECO[,c(1:8)]+1)),max(na.omit(subset_ELI_no_N_ECO[,c(1:8)]+1)))

clust.euclid.average <- hclust(dist(t(log2(subset_ELI_no_N_ECO[,c(1:8)]+1))),method="average")
#clust.euclid.average2 <- hclust(dist(t(na.omit(subset_ELI_NORM_BASIC_ECO[,c(1:8)]))),method="average")


pdf(file="QC_PME12_ECO_May19.pdf", width = 15, height = 15)
boxplot(log2(subset_ELI_no_N_ECO[,c(1:8)]+1),names=colnames(subset_ELI_no_N_ECO)[c(1:8)], which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(subset_ELI_no_N_ECO[,c(1:8)]+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
#boxplot(log2(subset_ELI_NORM_BASIC_ECO[,c(1:8)]+1), cex.axis=0.7, las=2)
#multi("density", subset_ELI_NORM_BASIC_ECO[,c(1:8)], xlim, "Density","","")
#plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_no_N_ECO[,c(1:8)], tipo, "DetOut_PME12_ECO_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_no_N_ECO[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_no_N_ECO[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_ECO_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()


acol = sample(brewer.pal(8, "Dark2"), ncol(subset_ELI_no_N_YEAST)-1, replace = (8 < ncol(subset_ELI_no_N_YEAST)-1))
xlim_r = c(min(na.omit(log2(subset_ELI_no_N_YEAST[,c(1:8)]+1))),max(na.omit(log2(subset_ELI_no_N_YEAST[,c(1:8)]+1))))
xlim = c(min(na.omit(subset_ELI_no_N_YEAST[,c(1:8)]+1)),max(na.omit(subset_ELI_no_N_YEAST[,c(1:8)]+1)))

clust.euclid.average <- hclust(dist(t(log2(subset_ELI_no_N_YEAST[,c(1:8)]+1))),method="average")
#clust.euclid.average2 <- hclust(dist(t(na.omit(subset_ELI_NORM_BASIC_YEAST[,c(1:8)]))),method="average")


pdf(file="QC_PME12_YEAST_May19.pdf", width = 15, height = 15)
boxplot(log2(subset_ELI_no_N_YEAST[,c(1:8)]+1),names=colnames(subset_ELI_no_N_YEAST)[c(1:8)], which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(subset_ELI_no_N_YEAST[,c(1:8)]+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
#boxplot(log2(subset_ELI_NORM_BASIC_YEAST[,c(1:8)]+1), cex.axis=0.7, las=2)
#multi("density", subset_ELI_NORM_BASIC_YEAST[,c(1:8)], xlim, "Density","","")
#plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_YEAST)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_no_N_YEAST[,c(1:8)], tipo, "DetOut_PME12_YEAST_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_no_N_YEAST[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_no_N_YEAST[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_YEAST_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()




















@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2






clust.euclid.average <- hclust(dist(t(log2(subset_ELI_no_N_BASIC_YEAST[,c(1:8)]+1))),method="average")
clust.euclid.average2 <- hclust(dist(t(na.omit(subset_ELI_NORM_BASIC_YEAST[,c(1:8)]))),method="average")


pdf(file="QC_PME12_YEAST_May19.pdf", width = 15, height = 15)
boxplot(log2(subset_ELI_no_N_BASIC_YEAST[,c(1:8)]+1),names=colnames(subset_ELI_no_N_BASIC_YEAST)[c(1:8)], which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(subset_ELI_no_N_BASIC_YEAST[,c(1:8)]+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(log2(subset_ELI_NORM_BASIC_YEAST[,c(1:8)]+1), cex.axis=0.7, las=2)
multi("density", subset_ELI_NORM_BASIC_YEAST[,c(1:8)], xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_YEAST)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_NORM_BASIC_YEAST[,c(1:8)], tipo, "DetOut_PME12_YEAST_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_NORM_BASIC_YEAST[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_NORM_BASIC_YEAST[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_YEAST_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

clust.euclid.average <- hclust(dist(t(log2(subset_ELI_no_N_BASIC_ECO[,c(1:8)]+1))),method="average")
clust.euclid.average2 <- hclust(dist(t(na.omit(subset_ELI_NORM_BASIC_ECO[,c(1:8)]))),method="average")



pdf(file="QC_PME12_ECO_May19.pdf", width = 15, height = 15)
boxplot(log2(subset_ELI_no_N_BASIC_ECO[,c(1:8)]+1),names=colnames(subset_ELI_no_N_BASIC_ECO)[c(1:8)], which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(subset_ELI_no_N_BASIC_ECO[,c(1:8)]+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(log2(subset_ELI_NORM_BASIC_ECO[,c(1:8)]+1), cex.axis=0.7, las=2)
multi("density", subset_ELI_NORM_BASIC_ECO[,c(1:8)], xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()

#tipo <- factor(unlist(sapply(colnames(subset_ELI_NORM_BASIC_ECO)[,c(1:8)], FUN=function(x) paste(unlist(strsplit(x, "_"))[c(1:8)],collapse="_"))))
tipo<-c("A","A","A","A","B","B","B","B")
a <- detOutliers(subset_ELI_NORM_BASIC_ECO[,c(1:8)], tipo, "DetOut_PME12_ECO_May19.pdf", 20, 20)
x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(subset_ELI_NORM_BASIC_ECO[,c(1:8)]))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo), size = 5) + geom_text(aes(label=colnames(subset_ELI_NORM_BASIC_ECO[,c(1:8)])), hjust=0, vjust=0, nudge_x=-5, nudge_y=3) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "PCA_PME12_ECO_Abr19.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()
