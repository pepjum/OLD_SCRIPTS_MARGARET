
name_sample_files<-list.files("/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/", pattern="dataPSMMat_filterPSM_id.txt")
file_path<-paste0("/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/",name_sample_files)


for (i in 1:length(file_path)){
    data<-read.csv2(file_path[i], header=F)
    if(i==1){
    results<-data.frame(data)
    }else{
        results<-rbind(results, data)
    }
}


results<-unique(results)

write.table(results, file="/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/id_list.txt", col.names=F, row.names=F, quote=F, sep="\t")
#GO TO swissPROT AND DOWNLOAD THE CORRELATION TABLE for IDs in NEXTPROT. SAVE AS id_table.txt



#changing swissprot id to Nextprot_id
mierda<-read.table("/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/id_table.txt", header=T, stringsAsFactors=F)
names(mierda)<-c("id", "NextProt")

name_sample_files2<-list.files("/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/", pattern="dataPSMMat_filterPSMNSAF_SC.rda")
file_path2<-paste0("/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/",name_sample_files2)

for (i in 1:length(file_path2)){
    cat(file_path2[i],"\n")
    data<-get(load(file_path2[i]))
    if(i==1){
    results<-data.frame(data)
    }else{
        results<-rbind(results, data)
    }
}

results$id<-sapply(strsplit(paste(results$ProteinAccession),"\\|"), "[", 2)

results_NX<-merge(results, mierda, by="id", all.x=T)

write.table(unique(paste(results_NX$NextProt)), file="/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/nextprot_id.txt", col.names=F, row.names=F, quote=F, sep="\t")

#ensembl<-read.table("/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/mart_export.txt", header=T)

nexprot_ensg<-read.table("/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/nextprot_ensg.txt", header=F)

results_NX_ensg<-merge(results_NX, nexprot_ensg, by.x="NextProt", by.y="V1", all.x=T)

final_results<-results_NX_ensg[,c("V2","sample","NSAF")]

matrix_ENSG_names<-data.frame("ENSG"=unique(paste(final_results$V2)))

names_final_results<-unique(paste(final_results$sample))


####hacerlo con matriz

the_matrix<-matrix(NA,nrow=length(unique(paste(final_results$V2))), ncol=length(names_final_results))
rownames(the_matrix)<-unique(paste(final_results$V2))
colnames(the_matrix)<-unique(paste(final_results$sample))

for(i in 1:nrow(final_results)){
    cat(i,"\n")
    the_matrix[paste(final_results$V2[i]),paste(final_results$sample[i])]<-as.numeric(final_results$NSAF[i])
}


saveRDS(the_matrix, file="/home/margaret/data/pepe/15_NCI60_uPE_ABR18/Results/matrix_ENSG_NCI60_quantification3.rds")
