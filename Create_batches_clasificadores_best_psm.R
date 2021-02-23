#### dividir en trozos
#unique(paste(tandem$PSM=(874060)  # esto deberian de ser las querys
#1522522/50 =30450

#nrow(30038694)
#mascot DP unique$query==(763601)
z<-1
i=1
n_cores<-31

mascot_querys<-length(unique(paste(results_mascot$Query)))
index_number<-mascot_querys/n_cores
index_number<-as.numeric(lapply(strsplit(paste(index_number),"\\."),"[",1))-1
mascot_folder<-"/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/Dat_Files/"

create_batches<-function(results_mascot,index_number,n_cores,search_engine, mascot_folder){
    for (j in seq(1:n_cores)){
        cat(j,"\n")
        if(z <= n_cores-1){
            z<-z+1
            kk<-unique(paste(results_mascot$Query))[i:(i+index_number)]
            cat(head(kk),"\n")
            batch<-results_mascot[which(paste(results_mascot$Query) %in% kk),]
            save(batch, file=paste0(mascot_folder,search_engine,"_batch_",j,".Rdata"))
            i<-i+index_number+1

        }else{
            i<-i+index_number+1
            kk<-unique(paste(results_mascot$Query))[i:length(unique(paste(results_mascot$Query)))]
            batch<-results_mascot[which(paste(results_mascot$Query) %in% kk),]
            save(batch, file=paste0(mascot_folder,search_engine,"batch_",j,".Rdata"))
        }
    }
}

batches_mascot<-create_batches(results_mascot, index_number,n_cores, "mascot", mascot_folder)

comet_querys<-length(unique(paste(results_comet$Query)))
index_number<-comet_querys/n_cores
index_number<-as.numeric(lapply(strsplit(paste(index_number),"\\."),"[",1))-1
comet_folder<-"/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/Comet_files/"

batches_comet<-create_batches(results_comet, index_number,n_cores, "comet", comet_folder)

omssa_querys<-length(unique(paste(results_omssa$Query)))
index_number<-omssa_querys/n_cores
index_number<-as.numeric(lapply(strsplit(paste(index_number),"\\."),"[",1))-1
omssa_folder<-"/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/Omssa_files/"

batches_omssa<-create_batches(results_omssa, index_number,n_cores, "omssa", omssa_folder)

tandem_querys<-length(unique(paste(results_tandem$Query)))
index_number<-tandem_querys/n_cores
index_number<-as.numeric(lapply(strsplit(paste(index_number),"\\."),"[",1))-1
tandem_folder<-"/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/Tandem_Files/"

batches_tandem<-create_batches(results_tandem, index_number,n_cores, "tandem", tandem_folder)


z<-1
i=1

for (j in seq(1:50)){
    cat(j,"\n")
    if(z <= 49){
        z<-z+1
        kk<-unique(paste(results_mascot$Query))[i:(i+17481)]
        cat(head(kk),"\n")
        batch<-results_mascot[which(paste(results_mascot$Query) %in% kk),]
        save(batch, file=paste0("tandem_batch_",j,".Rdata"))
        i<-i+17481+1

    }else{
        i<-i+17481+1
        kk<-unique(paste(results_mascot$Query))[i:length(unique(paste(results_mascot$Query)))]
        batch<-results_mascot[which(paste(results_mascot$Query) %in% kk),]
        save(batch, file=paste0("tandem_batch_",j,".Rdata"))
    }
}

##### lo mismo para comet

z<-1
i=1

for (j in seq(1:60)){
    cat(j,"\n")
    if(z <= 59){
        z<-z+1
        kk<-unique(paste(comet_file_load$Query))[i:(i+23594)]
        cat(head(kk),"\n")
        batch<-comet_file_load[which(paste(comet_file_load$Query) %in% kk),]
        save(batch, file=paste0("comet_batch_",j,".Rdata"))
        i<-i+23594+1

    }else{
        i<-i+23594+1
        kk<-unique(paste(comet_file_load$Query))[i:length(unique(paste(comet_file_load$Query)))]
        batch<-comet_file_load[which(paste(comet_file_load$Query) %in% kk),]
        save(batch, file=paste0("comet_batch_",j,".Rdata"))
    }
}

files_batch<-list.files(getwd(), pattern="comet_batch")
files_batched<-paste(getwd(),files_batch, sep="/")
write.table(files_batched,"files_to_best_psm_comet.txt", col.names=F,row.names=F, quote=F)


### batches matriz para PeptideMASS


for (j in seq(1:58)){
    cat(j,"\n")
    if(z <= 57){
        z<-z+1
        kk<-unique(paste(ALL_EXPERIMENT$PSM))[i:(i+38995)]
        cat(head(kk),"\n")
        batch<-ALL_EXPERIMENT[which(paste(ALL_EXPERIMENT$PSM) %in% kk),]
        save(batch, file=paste0("ALL_EXPERIMENT_batch_",j,".Rdata"))
        i<-i+38995+1

    }else{
        i<-i+38995+1
        kk<-unique(paste(ALL_EXPERIMENT$PSM))[i:length(unique(paste(ALL_EXPERIMENT$PSM)))]
        batch<-ALL_EXPERIMENT[which(paste(ALL_EXPERIMENT$PSM) %in% kk),]
        save(batch, file=paste0("ALL_EXPERIMENT_batch_",j,".Rdata"))
    }
}

for (j in seq(1:58)){
    cat(j,"\n")
    if(z <= 57){
        z<-z+1
        kk<-unique(paste(TMP_OBJECT$PSM))[i:(i+38995)]
        cat(head(kk),"\n")
        batch<-TMP_OBJECT[which(paste(TMP_OBJECT$PSM) %in% kk),]
        save(batch, file=paste0("TMP_OBJECT_batch_",j,".Rdata"))
        i<-i+38995+1

    }else{
        i<-i+38995+1
        kk<-unique(paste(TMP_OBJECT$PSM))[i:length(unique(paste(TMP_OBJECT$PSM)))]
        batch<-TMP_OBJECT[which(paste(TMP_OBJECT$PSM) %in% kk),]
        save(batch, file=paste0("TMP_OBJECT_batch_",j,".Rdata"))
    }
}

files_batch<-list.files(getwd(), pattern="ALL_EXPERIMENT_batch_")
files_batched<-paste(getwd(),files_batch, sep="/")
write.table(files_batched,"files_to_ALL_EXPERIMENT_PeptideMass.txt", col.names=F,row.names=F, quote=F)
