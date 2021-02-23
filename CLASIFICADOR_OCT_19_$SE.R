

library(protr)
library(caret)
library(ggplot2)    # features dependientes de la secuencia del peptido
library(gridExtra)
library(stringr)
library(plyr)
library(dplyr)

source("/home/margaret/data/pepe/scripts/functions_Classifier.R")

path<-args[1] #"/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/"

mascot_folder<-"Dat_Files/"
Comet_folder<-"Comet_files/"
Tandem_folder<-"Tandem_Files/"
Omssa_folder<-"Omssa_files/"



MGF_folder<-"MGFFiles"

mascot_folder<-paste0(path,mascot_folder)
Comet_folder<-paste0(path,Comet_folder)
Tandem_folder<-paste0(path,Tandem_folder)
MGF_folder<-paste0(path, MGF_folder)
Omssa_folder<-paste0(path, Omssa_folder)
#decoy_id<-args[5]    #"DECOY"
#minn_AA <- args[6]   # 9
#maxn_AA <- args[7]  # 30

mascot_files<-list.files(mascot_folder,pattern="dataPSMMat.rda")
Comet_files<-list.files(Comet_folder,pattern="dataPSMMat.rda")
Tandem_files<-list.files(Tandem_folder,pattern="dataPSMMat.rda")
Omssa_files<-list.files(Omssa_folder,pattern="dataPSMMat.rda")

results_mascot<-data.frame()
for(i in 1:length(mascot_files)){
    cat("loading ",mascot_files[i],"\n")
    mascot_file_load<-get(load(paste0(mascot_folder,mascot_files[i])))
    results_mascot<-rbind(results_mascot,mascot_file_load)
}

### ME quedo solo con los best psm
mascot_querys<-length(unique(paste(results_mascot$Query)))
index_number<-mascot_querys/n_cores
index_number<-as.numeric(lapply(strsplit(paste(index_number),"\\."),"[",1))-1


n_cores<-55

create_batches<-function(results_mascot,index_number,n_cores,search_engine, mascot_folder){
    z<-1
    i=1

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
            save(batch, file=paste0(mascot_folder,search_engine,"_batch_",j,".Rdata"))
        }
    }
}

batches_mascot<-create_batches(results_mascot, index_number,n_cores, "mascot", mascot_folder)



results_comet<-data.frame()
for(i in 1:length(Comet_files)){
    cat("loading ",Comet_files[i],"\n")
    comet_file_load<-get(load(paste0(Comet_folder,Comet_files[i])))
    results_comet<-rbind(results_comet,comet_file_load)
}

comet_querys<-length(unique(paste(results_comet$Query)))
index_number<-comet_querys/n_cores
index_number<-as.numeric(lapply(strsplit(paste(index_number),"\\."),"[",1))-1

batches_comet<-create_batches(results_comet, index_number,n_cores, "comet", Comet_folder)

results_Omssa<-data.frame()
for(i in 1:length(Omssa_files)){
    cat("loading ",Omssa_files[i],"\n")
    Omssa_file_load<-get(load(paste0(Omssa_folder,Omssa_files[i])))
    cat(names(Omssa_file_load),"\n")
    results_Omssa<-rbind(results_Omssa,Omssa_file_load)
}

omssa_querys<-length(unique(paste(results_Omssa$Query)))
index_number<-omssa_querys/n_cores
index_number<-as.numeric(lapply(strsplit(paste(index_number),"\\."),"[",1))-1


batches_Omssa<-create_batches(results_Omssa, index_number,n_cores, "Omssa", Omssa_folder)

tandem_querys<-length(unique(paste(results_tandem$Query)))
index_number<-tandem_querys/n_cores
index_number<-as.numeric(lapply(strsplit(paste(index_number),"\\."),"[",1))-1



results_tandem<-data.frame()
for(i in 1:length(Tandem_files)){
    cat("loading ",Tandem_files[i],"\n")
    tandem_file_load<-get(load(paste0(Tandem_folder,Tandem_files[i])))
    cat(names(tandem_file_load),"\n")
    results_tandem<-rbind(results_tandem,tandem_file_load)
}

batches_Tandem<-create_batches(results_tandem, index_number,n_cores, "Tandem", Tandem_folder)


############### LANZAR EL LAUNCH_BEST_PSM #para cada uno de los motores de busqueda

y<-paste('sh /home/margaret/data/pepe/scripts/launch_best_psm.sh -f ',paste(Comet_folder))
system(y)



# mascot_best_psm<-data.frame()
# for (query in unique(paste(mascot_file_load$Query))){
#     print(query)
#     tmp<-mascot_file_load[which(mascot_file_load$Query==query),]
#     if((nrow(tmp) >=2 )& length(unique(paste(tmp$score) >1))){
#         tmp_selected<-tmp[which.max(as.numeric(paste(tmp$score))),]
#         mascot_best_psm<-rbind(mascot_best_psm,tmp_selected)
#     }else if(nrow(tmp)==1){
#         mascot_best_psm<-rbind(mascot_best_psm,tmp)
#     }
#
# }
