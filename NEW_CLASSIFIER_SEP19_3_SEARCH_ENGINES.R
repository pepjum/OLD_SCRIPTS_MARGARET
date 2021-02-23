args=(commandArgs(TRUE))

library(protr)
library(caret)
library(ggplot2)    # features dependientes de la secuencia del peptido
library(gridExtra)
library(stringr)
library(plyr)
library(dplyr)

source("/home/margaret/data/pepe/scripts/functions_Classifier.R")

path<-args[1] #"/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/"

mascot_folder<-"Dat_Files"
Comet_folder<-"Comet_files"
Tandem_folder<-"Tandem_Files"
Omssa_folder<-"Omssa_files"

MGF_folder<-"MGFFiles"

mascot_folder<-paste0(path,mascot_folder)
Comet_folder<-paste0(path,Comet_folder)
Tandem_folder<-paste0(path,Tandem_folder)
MGF_folder<-paste0(path, MGF_folder)
#decoy_id<-args[5]    #"DECOY"
#minn_AA <- args[6]   # 9
#maxn_AA <- args[7]  # 30

mascot_file<-list.files(mascot_folder,pattern="dataPSMMat.rda")
Comet_file<-list.files(Comet_folder,pattern="dataPSMMat.rda")
Tandem_file<-list.files(Tandem_folder,pattern="dataPSMMat.rda")
Omssa_file<-list.files(Omssa_folder,pattern="dataPSMMat.rda")

mascot_file_load<-get(load(paste(mascot_folder,mascot_file,sep="/")))
mascot_file_load$Scan_tmp<-lapply(strsplit(paste(mascot_file_load$ScanString)," "),"[",1)
mascot_file_load$Scan_tmp2<-lapply(strsplit(paste(mascot_file_load$Scan_tmp),"\\."),"[",2)
mascot_file_load$Query_tmp<-lapply(strsplit(paste(mascot_file_load$Scan_tmp),"\\."),"[",1)
mascot_file_load$ScanNumber<-mascot_file_load$Scan_tmp2
mascot_file_load$Query<-paste(mascot_file_load$Query_tmp,mascot_file_load$Scan_tmp2,sep="_")
mascot_file_load$Scan_tmp<-NULL
mascot_file_load$Scan_tmp2<-NULL
mascot_file_load$Query_tmp<-NULL
mascot_file_load$Rank<-NULL
mascot_file_load$PSM<-paste(mascot_file_load$Query, mascot_file_load$PeptideSeq, sep="_")
mascot_file_load<-unique(mascot_file_load)
mascot_file_load$PSM<-str_replace(paste(mascot_file_load$PSM), "%2d", "-")

#primero me quedo con el mejor PSM de cada espectro
mascot_best_psm<-data.frame()
for (query in unique(paste(mascot_file_load$Query))){
    print(query)
    tmp<-mascot_file_load[which(mascot_file_load$Query==query),]
    if((nrow(tmp) >=2 )& length(unique(paste(tmp$score) >1))){
        tmp_selected<-tmp[which.max(as.numeric(paste(tmp$score))),]
        mascot_best_psm<-rbind(mascot_best_psm,tmp_selected)
    }else if(nrow(tmp)==1){
        mascot_best_psm<-rbind(mascot_best_psm,tmp)
    }

}
#### OTRA OPCION más rapida, quitar el ProteinAccession, unique, y despues PeptideMatch cuando haya que hacer el ProtFDR
mascot_dummy_code<-mascot_file_load[,c("PeptideSeq","ProteinAccession")]
save(mascot_dummy_code, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/mascot_Protein_Accessions.Rdata")

mascot_best_psm$ProteinAccession<-NULL
mascot_best_psm<-unique(mascot_best_psm)
##### OJO! el SCAN de COMET son es el index del mgf, no es el scanNumber. Cambiar del objeto de entrada, y recuperar del mgf el SCAN antes de seguir
####   PREPROCESAR COMET. ESTO DESPUES MEJOR SACARLO A OTRO SCRIPT)
save(mascot_best_psm, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/results_mascot_best_psm.Rdata")

comet_file_load<-get(load(paste(Comet_folder, Comet_file, sep="/")))
comet_file_load$Fraction<-lapply(strsplit(paste(comet_file_load$datfile),"_"),"[",1)
comet_file_load$MGFfile<-paste0(comet_file_load$Fraction,".mgf")
comet_file_load$Query_tmp<-NULL
comet_file_load$mgfIndex<-comet_file_load$ScanNumber
comet_file_load$ScanNumber<-NA
comet_file_load$dummyCode<-paste(comet_file_load$Fraction,comet_file_load$mgfIndex,sep="_")


####
samples<-list.dirs(MGF_folder)
samples<-samples[-1]

trazability_out<-data.frame()

for(i in 1:length(samples)){

    files<-list.files(samples[i], pattern=".mgf")
    files<-paste(samples[i], files, sep="/")

    for (j in 1:length(files)){
        file_out<-paste0(lapply(strsplit(files[j],"\\."),"[",1),".index")

        y<- paste0('grep "^TITLE" ', files[j],' > ', file_out)
        cat(y,"\n")
        system(y)
    }

    files_index_comet<-list.files(samples[i], pattern=".index")
    files_index_comet<-paste(samples[i],files_index_comet, sep="/")


    for(j in 1:length(files_index_comet)){
        cat(files_index_comet[j],"\n")
        trazability<-read.table(files_index_comet[j], header=F)
        fraction_tmp<-paste(lapply(strsplit(paste(files_index_comet[j]),"\\/"),"[", 9))
        fraction_tmp2<-sapply(strsplit(paste(fraction_tmp),"\\."),"[",1)
        df_tmp<-data.frame("scanstring"=trazability[,5])
        names(df_tmp)<-"scanstring"
        df_tmp$mgfIndexes<-seq(1:nrow(df_tmp))
        df_tmp$fraction_t<-fraction_tmp2
        df_tmp$dummyCode<-paste(paste(as.character(df_tmp$fraction_t)),paste(as.character(df_tmp$mgfIndexes)),sep="_")
        df_tmp<-df_tmp[,c(1,4)]
        trazability_out<-rbind(trazability_out, df_tmp)
    }
}


#### merge objeto trazabilidad, con el de comet

Comet_merged<-merge(comet_file_load,trazability_out, by="dummyCode", all.x=TRUE)
Comet_merged<-unique(Comet_merged)
Comet_merged$ScanNumber<-gsub("[^0-9.]", "",  Comet_merged$scanstring)
Comet_merged$PSM<-paste(Comet_merged$Fraction, Comet_merged$ScanNumber, Comet_merged$PeptideSeq,sep="_")
Comet_merged$Query<-paste(Comet_merged$Fraction, Comet_merged$ScanNumber, sep="_")

save(Comet_merged, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Comet_files/results_Peptides_comet_PXD001383_dataPSMMat_corrected.rda")

#quedarme con el mejor PSM para cada espectro

###### ir al script de Create_batches_clasificadores.R para dividir y hacerlo mas rapido. Una vez hecho, ejecutar best_PSM.R para cada batch

# comet_best_psm<-data.frame()
# for (query in unique(paste(Comet_merged$Query))){
#     print(query)
#     tmp<-Comet_merged[which(Comet_merged$Query==query),]
#     if((nrow(tmp) >=2 )& length(unique(paste(tmp$score) >1))){
#         tmp_selected<-tmp[which.max(as.numeric(paste(tmp$score))),]
#         comet_best_psm<-rbind(comet_best_psm,tmp_selected)
#     }else if(nrow(tmp)==1){
#         comet_best_psm<-rbind(comet_best_psm,tmp)
#     }
#
# }




######  elimino duplicidades de cada PSM por haber asignado proteinas distintas para el proteingroupping
comet_dummy_code<-Comet_merged[,c("PeptideSeq","ProteinAccession")]
save(comet_dummy_code, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/comet_Protein_Accessions.Rdata")

comet_best_psm$ProteinAccession<-NULL
comet_best_psm<-unique(comet_best_psm)

##### load Tandem
tandem_file_load<-get(load(paste(Tandem_folder,Tandem_file,sep="/")))
tandem_file_load$Query_tmp<-lapply(strsplit(paste(tandem_file_load$ScanString),"\\."),"[",1)
tandem_file_load$Query<-paste(tandem_file_load$Query_tmp,tandem_file_load$ScanNumber,sep="_")
tandem_file_load$Query_tmp<-NULL
tandem_file_load$PSM<-paste(tandem_file_load$Query, tandem_file_load$PeptideSeq, sep="_")

#primero me quedo con el mejor PSM de cada espectro en TANDEM


tandem_best_psm<-data.frame()
for (query in unique(paste(tandem_file_load$Query))){
    print(query)
    tmp<-tandem_file_load[which(tandem_file_load$Query==query),]
    if((nrow(tmp) >=2 )& length(unique(paste(tmp$score) >1))){
        tmp_selected<-tmp[which.max(as.numeric(paste(tmp$score))),]
        tandem_best_psm<-rbind(tandem_best_psm,tmp_selected)
    }else if(nrow(tmp)==1){
        tandem_best_psm<-rbind(tandem_best_psm,tmp)
    }

}


tandem_dummy_code<-tandem_file_load[,c("PeptideSeq","ProteinAccession")]
save(tandem_dummy_code, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/tandem_Protein_Accessions.Rdata")


tandem_best_psm$ProteinAccession<-NULL
tandem_best_psm<-unique(tandem_best_psm)

##### CREAR EL DATAFRAME MATRIZ

unique_mascot_PSM<-unique(paste(mascot_best_psm$PSM))
unique_comet_PSM<-unique(paste(comet_best_psm$PSM))
unique_tandem_PSM<-unique(paste(tandem_best_psm$PSM))

list_querys<-c(unique_mascot_PSM,unique_comet_PSM, unique_tandem_PSM)
list_querys_unique<-unique(list_querys)  #2261740

mascot_tmp<-mascot_best_psm[,c(18,4,6,8,13)]
names(mascot_tmp)<-c("PSM","PeptideMass_mascot","PeptideSeq_mascot","score_mascot","database_mascot")
mascot_tmp<-unique(mascot_tmp)
#mascot_tmp$psmFDR_mascot<-psmFDR(mascot_tmp, "score_mascot", "DECOY", 0, "database_mascot")

comet_tmp<-comet_best_psm[,c(23,6,8,10,15)]
names(comet_tmp)<-c("PSM","PeptideMass_comet","PeptideSeq_comet","score_comet","database_comet")
comet_tmp<-unique(comet_tmp)
#comet_tmp$psmFDR_comet<-psmFDR(comet_tmp, "score_comet", "DECOY", 0, "database_comet")

tandem_tmp<-tandem_best_psm[,c(20,5,7,14,15)]
names(tandem_tmp)<-c("PSM","PeptideMass_tandem","PeptideSeq_tandem","score_tandem","database_tandem")
tandem_tmp<-unique(tandem_tmp)
#tandem_tmp$psmFDR_tandem<-psmFDR(comet_tmp, "score_tandem", "DECOY", 0, "database_tandem")
rm(mascot_file_load, Comet_file_load, Comet_merged, tandem_file_load)



Query_df<-data.frame("PSM"=list_querys_unique)

Query_df<-merge(Query_df, mascot_tmp, by.x="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, comet_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, tandem_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, omssa_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)

Query_df$PeptideSeq<-lapply(strsplit(paste(Query_df$PSM),"_"),"[", 8)
rm(list_querys, list_querys_unique)


save(Query_df, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/MATRIX_PSM_initial.Rdata")
### ordenar por cada score, y hacer un sumatorio de todas las posiciones, para rankear por el sumatorio de las tres columnas
#### psmFDR para cada motor

cat("computing some features...\n")

# Query_df$PepLen_mascot<-as.numeric(paste(nchar(as.character(Query_df$PeptideSeq_mascot))))
# Query_df$PepLen_comet<-as.numeric(paste(nchar(as.character(Query_df$PeptideSeq_comet))))
# Query_df$PepLen_tandem<-as.numeric(paste(nchar(as.character(Query_df$PeptideSeq_tandem))))


Query_df$database_mascot<-as.character(paste(Query_df$database_mascot))
Query_df$database_comet<-as.character(paste(Query_df$database_comet))
Query_df$database_tandem<-as.character(paste(Query_df$database_tandem))

for(number in c(2,4,6,8,10,12)){
    cat(number,"\n")
    Query_df[,number]<-as.numeric(paste(Query_df[,number]))
}

cat("preparing training set....\n")

ALL_DECOY<-rbind(Query_df[which(Query_df$database_mascot=="D"),],Query_df[which(Query_df$database_comet=="D"),],Query_df[which(Query_df$database_tandem=="D"),])
# ALL_DECOY$database<-"DECOY"
# ALL_DECOY$database_mascot<-NULL
# ALL_DECOY$database_comet<-NULL
# ALL_DECOY$database_tandem<-NULL
ALL_DECOY<-unique(ALL_DECOY)

save(ALL_DECOY, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/ALL_DECOY.Rdata")

DECOYS_for_TRAINING<-sample_n(ALL_DECOY, nrow(ALL_DECOY)/2)

ALL_TARGET<-rbind(Query_df[which(Query_df$database_mascot=="T"),], Query_df[which(Query_df$database_comet=="T"),], Query_df[which(Query_df$database_tandem=="T"),])
# ALL_TARGET$database<-"TARGET"
# ALL_TARGET$database_mascot<-NULL
# ALL_TARGET$database_comet<-NULL
# ALL_TARGET$database_tandem<-NULL
ALL_TARGET<-unique(ALL_TARGET)

save(ALL_TARGET, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/ALL_TARGET.Rdata")

rm(Query_df)

ALL_EXPERIMENT<-rbind(ALL_TARGET,ALL_DECOY)    #### OBJETO CON EL EXPERIMENTO YA PREPARADO
ALL_EXPERIMENT<-unique(ALL_EXPERIMENT) #nrow(ALL_EXPERIMENT)

TMP_OBJECT<-rbind(ALL_TARGET, DECOYS_for_TRAINING)    #ESTE objecto es sobre el que hay que hacer primer FDR
TMP_OBJECT<-unique(TMP_OBJECT)    #nrow(TMP_OBJECT)

### subset de este objeto para cada motor para poder calcular el FDR

TMP_OBJECT_mascot_selected<-TMP_OBJECT[which(!is.na(TMP_OBJECT$score_mascot)),]
mascot_DF_with_FDR<-psmFDR(TMP_OBJECT_mascot_selected, "score_mascot", "D", concat_decoy=0, prot_id="database_mascot")
mascot_TMP_FDR<-mascot_DF_with_FDR[,c("PSM","psmFDR")]
names(mascot_TMP_FDR)[2]<-"psmFDR_mascot"

TMP_OBJECT<-left_join(TMP_OBJECT, mascot_TMP_FDR, by = "PSM")

# TMP_OBJECT<-merge(TMP_OBJECT, mascot_TMP_FDR, by="PSM", all.x=TRUE)
TMP_OBJECT<-unique(TMP_OBJECT)

TMP_OBJECT_comet_selected<-TMP_OBJECT[which(!is.na(TMP_OBJECT$score_comet)),]
comet_DF_with_FDR<-psmFDR(TMP_OBJECT_comet_selected, "score_comet", "D", concat_decoy=0, prot_id="database_comet")
comet_TMP_FDR<-comet_DF_with_FDR[,c("PSM","psmFDR")]
names(comet_TMP_FDR)[2]<-"psmFDR_comet"

TMP_OBJECT<-merge(TMP_OBJECT, comet_TMP_FDR, by="PSM", all.x=TRUE)

#TMP_OBJECT<-full_join(TMP_OBJECT, comet_TMP_FDR, by = "PSM")

TMP_OBJECT<-unique(TMP_OBJECT)

TMP_OBJECT_tandem_selected<-TMP_OBJECT[which(!is.na(TMP_OBJECT$score_tandem)),]
tandem_DF_with_FDR<-psmFDR(TMP_OBJECT_tandem_selected, "score_tandem", "D", concat_decoy=0, prot_id="database_tandem")
tandem_TMP_FDR<-tandem_DF_with_FDR[,c("PSM","psmFDR")]
names(tandem_TMP_FDR)[2]<-"psmFDR_tandem"

TMP_OBJECT<-merge(TMP_OBJECT, tandem_TMP_FDR, by="PSM", all.x=TRUE)
TMP_OBJECT<-unique(TMP_OBJECT)

save(TMP_OBJECT, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/TMP_OBJECT_PSM.Rdata")   #objeto para hacer rerank

rm(mascot_TMP_FDR,mascot_DF_with_FDR,comet_DF_with_FDR,comet_TMP_FDR,tandem_DF_with_FDR,tandem_TMP_FDR)

DF_FDR_PASSED_mascot<-TMP_OBJECT[which(na.omit(paste(TMP_OBJECT$psmFDR_mascot))<=0.01),]
DF_FDR_PASSED_comet<-TMP_OBJECT[which(na.omit(paste(TMP_OBJECT$psmFDR_comet))<=0.01),]
DF_FDR_PASSED_tandem<-TMP_OBJECT[which(na.omit(paste(TMP_OBJECT$psmFDR_tandem))<=0.01),]

DF_FDR_PASSED<-unique(rbind(DF_FDR_PASSED_comet,DF_FDR_PASSED_mascot,DF_FDR_PASSED_tandem))

save(DF_FDR_PASSED, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/DF_FDR_PASSED_SET_FOR_FIRST_ITERATION.Rdata")   # con este objeto + mitad decoys se entrena en cada iteracion, rerankeandolo en cada iteracion


#### PREPARAR LOS DATOS ANTES DE ENTRAR AL CLASIFICADOR
if(nrow(DF_FDR_PASSED==0)){

    cat("There are some DECOY PSMs hits with scores higher than the best TARGET PSM. Could not create training dataset. Please remove these hits and try again...\n")
    break
}

TRAINING_OBJECT<-DF_FDR_PASSED
rm(DF_FDR_PASSED)

TRAINING_OBJECT$psmFDR_mascot<-NULL
TRAINING_OBJECT$psmFDR_comet<-NULL
TRAINING_OBJECT$psmFDR_tandem<-NULL

TRAINING_OBJECT<-rbind(TRAINING_OBJECT,DECOYS_for_TRAINING)
save(TRAINING_OBJECT,file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/TRAINING_SET_FOR_FIRST_ITERATION.Rdata")

#### TERMINAR DE ARREGLAR OBJETO
cat("preparing ALL Experiment matrix","\n")
ALL_EXPERIMENT$PeptideSeq<-lapply(strsplit(paste(ALL_EXPERIMENT$PSM),"_"),"[",3)
ALL_EXPERIMENT$dummy<-seq(1:nrow(ALL_EXPERIMENT))
ALL_EXPERIMENT$PeptideSeq_comet<-NULL
ALL_EXPERIMENT$PeptideSeq_tandem<-NULL
ALL_EXPERIMENT$PeptideSeq_mascot<-NULL

ncores<-58
index<-nrow(ALL_EXPERIMENT)/ncores
index<-as.numeric(strsplit(paste(index),"\\.")[[1]][1]) -1

i<-1
z<-1

for (j in seq(1:ncores)){
    cat(j,"\n")
    if(z <= (ncores-1)){
        z<-z+1
        batch<-ALL_EXPERIMENT[i:(i+index),]
        cat("rows",nrow(batch),"\n")
        save(batch, file=paste0(getwd(),"/","ALL_EXPERIMENT_TEMP/",paste0("ALL_EXPERIMENT_batch_",j,".Rdata")))
        i<-i+index+1
        cat("i ",i,"\n")

    }else{
        cat("i",i,"\n")
        batch<-ALL_EXPERIMENT[i:nrow(ALL_EXPERIMENT),]
        cat("rows",nrow(batch),"\n")
        save(batch, file=paste0(getwd(),"/","ALL_EXPERIMENT_TEMP/",paste0("ALL_EXPERIMENT_batch_",j,".Rdata")))
    }
}

files_batch<-list.files(paste0(getwd(),"/ALL_EXPERIMENT_TEMP/"), pattern="ALL_EXPERIMENT_batch_")
files_batched<-paste(paste0(getwd(),"/ALL_EXPERIMENT_TEMP/"),files_batch, sep="/")
write.table(files_batched,"files_to_ALL_EXPERIMENT_PeptideMass.txt", col.names=F,row.names=F, quote=F)

y<-paste0('bash /home/margaret/data/pepe/scripts/launch_PeptideMass_and_Db_ALL_EXPERIMENT.sh -f ', paste0(getwd(),"/"))

system(y)

### load dummy objects

files_dummy<-list.files("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/ALL_EXPERIMENT_TEMP/", pattern="AND_DB.Rdata")
files_dummy<-paste0(getwd(),"/","ALL_EXPERIMENT_TEMP/",files_dummy)

dummy_all<-data.frame()
for (file in files_dummy){
    cat(file,"\n")
    dummy_load<-get(load(file))
    dummy_all<-rbind(dummy_all,dummy_load)
}

ALL_EXPERIMENT<-merge(ALL_EXPERIMENT,dummy_all, by="dummy", all.x=T)

ALL_EXPERIMENT<-ALL_EXPERIMENT[,c("PSM","PeptideSeq","PeptideMass","score_mascot","score_comet","score_tandem","database")]
ALL_EXPERIMENT$score_mascot<-as.numeric(paste(ALL_EXPERIMENT$score_mascot))
ALL_EXPERIMENT$score_comet<-as.numeric(paste(ALL_EXPERIMENT$score_comet))
ALL_EXPERIMENT$score_tandem<-as.numeric(paste(ALL_EXPERIMENT$score_tandem))

save(ALL_EXPERIMENT,file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/ALL_EXPERIMENT_MATRIX_with_NA_PREPARED.Rdata")


ALL_EXPERIMENT[is.na(ALL_EXPERIMENT)]<-0

save(ALL_EXPERIMENT, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/ALL_EXPERIMENT_MATRIX_PREPARED.Rdata")



cat("preparing training object","\n")

TRAINING_OBJECT$PeptideSeq<-lapply(strsplit(paste(TRAINING_OBJECT$PSM),"_"),"[",3)
TRAINING_OBJECT$PeptideSeq_comet<-NULL
TRAINING_OBJECT$PeptideSeq_mascot<-NULL
TRAINING_OBJECT$PeptideSeq_tandem<-NULL
TRAINING_OBJECT<-unique(TRAINING_OBJECT)

TRAINING_OBJECT$dummy<-seq(1:nrow(TRAINING_OBJECT))

ncores<-58

index<-nrow(TRAINING_OBJECT)/ncores
index<-as.numeric(strsplit(paste(index),"\\.")[[1]][1]) -1

i<-1
z<-1

for (j in seq(1:ncores)){
    cat(j,"\n")
    if(z <= (ncores-1)){
        z<-z+1
        batch<-TRAINING_OBJECT[i:(i+index),]
        cat("rows",nrow(batch),"\n")
        save(batch, file=paste0(getwd(),"/","TRAINING_TEMP/",paste0("TRAINING_OBJECT_batch_",j,".Rdata")))
        i<-i+index+1
        cat("i ",i,"\n")

    }else{
        cat("i",i,"\n")
        batch<-TRAINING_OBJECT[i:nrow(TRAINING_OBJECT),]
        cat("rows",nrow(batch),"\n")
        save(batch, file=paste0(getwd(),"/","TRAINING_TEMP/",paste0("TRAINING_OBJECT_batch_",j,".Rdata")))
    }
}

files_batch<-list.files(paste0(getwd(),"/TRAINING_TEMP/"), pattern="TRAINING_OBJECT_batch_")
files_batched<-paste(paste0(getwd(),"/TRAINING_TEMP/"),files_batch, sep="/")
write.table(files_batched,"files_to_TRAINING_OBJECT_PeptideMass.txt", col.names=F,row.names=F, quote=F)

y<-paste0('bash /home/margaret/data/pepe/scripts/launch_PeptideMass_and_Db_TMP_OBJECT.sh -f ', paste0(getwd(),"/"))

system(y)
###### WAIT to the end

files_dummy<-list.files("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/TRAINING_TEMP/", pattern="AND_DB.Rdata")
files_dummy<-paste0(getwd(),"/","TRAINING_TEMP/",files_dummy)
dummy_all<-data.frame()
for (i in seq(1:length(files_dummy))){
    cat(i,"\n")
    if(i==1){
    dummy_load<-get(load(files_dummy[i]))
    dummy_all<-dummy_load
    }else{
        dummy_load<-get(load(files_dummy[i]))
        dummy_all<-rbind(dummy_all,dummy_load)
    }
}

TRAINING_OBJECT<-merge(TRAINING_OBJECT,dummy_all, by="dummy", all.x=T)
TRAINING_OBJECT<-TRAINING_OBJECT[,c("PSM","PeptideSeq","PeptideMass","score_mascot","score_comet","score_tandem","database")]

TRAINING_OBJECT$score_mascot<-as.numeric(paste(TRAINING_OBJECT$score_mascot))
TRAINING_OBJECT$score_comet<-as.numeric(paste(TRAINING_OBJECT$score_comet))
TRAINING_OBJECT$score_tandem<-as.numeric(paste(TRAINING_OBJECT$score_tandem))

save(TRAINING_OBJECT, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/TRAINING_OBJECT__with_NAs_PREPARED.Rdata")


TRAINING_OBJECT[is.na(TRAINING_OBJECT)]<-0

save(TRAINING_OBJECT, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/TRAINING_OBJECT_PREPARED.Rdata")


cat("preparing decoys for training object","\n")

DECOYS_for_TRAINING$PeptideSeq<-lapply(strsplit(paste(DECOYS_for_TRAINING$PSM),"_"),"[",3)
DECOYS_for_TRAINING$PeptideSeq_comet<-NULL
DECOYS_for_TRAINING$PeptideSeq_mascot<-NULL
DECOYS_for_TRAINING$PeptideSeq_tandem<-NULL
DECOYS_for_TRAINING<-unique(DECOYS_for_TRAINING)

DECOYS_for_TRAINING$dummy<-seq(1:nrow(DECOYS_for_TRAINING))

ncores<-58

index<-nrow(DECOYS_for_TRAINING)/ncores
index<-as.numeric(strsplit(paste(index),"\\.")[[1]][1]) -1

i<-1
z<-1

for (j in seq(1:ncores)){
    cat(j,"\n")
    if(z <= (ncores-1)){
        z<-z+1
        batch<-DECOYS_for_TRAINING[i:(i+index),]
        cat("rows",nrow(batch),"\n")
        save(batch, file=paste0(getwd(),"/","DECOYS_for_TRAINING_TEMP/",paste0("DECOYS_for_TRAINING_batch_",j,".Rdata")))
        i<-i+index+1
        cat("i ",i,"\n")

    }else{
        cat("i",i,"\n")
        batch<-DECOYS_for_TRAINING[i:nrow(DECOYS_for_TRAINING),]
        cat("rows",nrow(batch),"\n")
        save(batch, file=paste0(getwd(),"/","DECOYS_for_TRAINING_TEMP/",paste0("DECOYS_for_TRAINING_batch_",j,".Rdata")))
    }
}

files_batch<-list.files(paste0(getwd(),"/DECOYS_for_TRAINING_TEMP/"), pattern="DECOYS_for_TRAINING_batch_")
files_batched<-paste(paste0(getwd(),"/DECOYS_for_TRAINING_TEMP/"),files_batch, sep="/")
write.table(files_batched,"files_to_DECOYS_for_TRAINING_PeptideMass.txt", col.names=F,row.names=F, quote=F)

y<-paste0('bash /home/margaret/data/pepe/scripts/launch_PeptideMass_and_Db_DECOYS.sh -f ', paste0(getwd(),"/"))

system(y)

files_dummy<-list.files("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/DECOYS_for_TRAINING_TEMP/", pattern="AND_DB.Rdata")
files_dummy<-paste0(getwd(),"/","DECOYS_for_TRAINING_TEMP/",files_dummy)

dummy_all<-data.frame()
for (file in files_dummy){
    cat(file,"\n")
    dummy_load<-get(load(file))
    dummy_all<-rbind(dummy_all,dummy_load)
}

DECOYS_for_TRAINING<-merge(DECOYS_for_TRAINING,dummy_all, by="dummy", all.x=T)
DECOYS_for_TRAINING<-DECOYS_for_TRAINING[,c("PSM","PeptideSeq","PeptideMass","score_mascot","score_comet","score_tandem","database")]

DECOYS_for_TRAINING$score_mascot<-as.numeric(paste(DECOYS_for_TRAINING$score_mascot))
DECOYS_for_TRAINING$score_comet<-as.numeric(paste(DECOYS_for_TRAINING$score_comet))
DECOYS_for_TRAINING$score_tandem<-as.numeric(paste(DECOYS_for_TRAINING$score_tandem))

DECOYS_for_TRAINING[is.na(DECOYS_for_TRAINING)]<-0

save(DECOYS_for_TRAINING, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/DECOYS_for_TRAINING_PREPARED.Rdata")




set.seed(1)
ctrl <- trainControl(method="repeatedcv", repeats=3, classProbs=TRUE, summaryFunction= twoClassSummary,savePredictions = TRUE,verboseIter=TRUE)

iterations<-5
############################################ AQUI EMPIEZA EL CLASIFICADOR

for (i in seq(1:iteration_n)){

    if (i==1){
        cat("training method iteration ", i,"\n")
        #glmFit<-train(database~., data=TRAINING_OBJECT[,c(3:6)], method="glm", family=binomial, trControl=ctrl, metric="ROC")
        glmFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "glm", metric="ROC", trControl = ctrl)

        cat("reranking test data", i,"\n")

        iteration<-predict(glmFit, TMP_OBJECT[,c(3:6)], type = 'prob')  #primera
        TMP_OBJECT[,8]<-as.numeric(iteration[,2])
        names(TMP_OBJECT)[8]<-paste0("p_TARGET_",i)
        a<-names(TMP_OBJECT[20])
        FIRST_ITERATION_FDR<-psmFDR(TMP_OBJECT, a, "DECOY",concat_decoy=0, prot_id="Proteins")
        df_temp<-FIRST_ITERATION_FDR_passed[,c(1:19)]

        FIRST_ITERATION_FDR_passed<-FIRST_ITERATION_FDR[which(FIRST_ITERATION_FDR$psmFDR <=0.01),]
        dataframes[[counter]]<-df_temp
        psms[counter]<-nrow(FIRST_ITERATION_FDR_passed)
        peptides[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed$Peptide)))
        proteins[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed$Proteins)))
        glms[[counter]]<-glmFit

        cat("applying to full data ",i,"\n")

        iteration_all<-predict(glmFit, dataPSMMat[,c(5,6,7,10,19)], type = 'prob')  #primera
        df[,20]<-as.numeric(iteration_all[,2])
        names(df)[20]<-paste0("p_TARGET_",i)
        a<-names(df[20])
        FIRST_ITERATION_FDR_all<-psmFDR(df, a, "DECOY",concat_decoy=0, prot_id="database")
        FIRST_ITERATION_FDR_passed_all<-FIRST_ITERATION_FDR_all[which(FIRST_ITERATION_FDR_all$psmFDR <=0.01),]
        df_temp_all<-FIRST_ITERATION_FDR_passed_all[,c(1:19)]

        dataframes_all[[counter]]<-df_temp_all
        psms_all[counter]<-nrow(FIRST_ITERATION_FDR_passed_all)
        peptides_all[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed_all$Peptide)))
        proteins_all[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed_all$Proteins)))

        only_targets<-FIRST_ITERATION_FDR_passed_all[which(FIRST_ITERATION_FDR_passed_all$Label=="TARGET"),]
        psms_all_clean[counter]<-nrow(only_targets)
        peptides_all_clean[counter]<-length(unique(paste(only_targets$Peptide)))
        proteins_all_clean[counter]<-length(unique(paste(only_targets$Proteins)))

        glms_all[[counter]]<-glmFit

        counter<-counter+1

    }else{

        trainDescr_tmp<-rbind(dataframes[[counter-1]], DECOYS_for_TRAINING)

        cat("training method iteration ", i,"\n")

        glmFit_2<-train(database~., data=trainDescr_tmp[,c(5,6,7,10,19)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

        glms[[counter]]<-glmFit_2
        glms_all[[counter]]<-glmFit_2

        cat("reranking test data", i,"\n")

        iteration<-predict(glmFit_2, TMP_OBJECT[,c(5,6,7,10,19)], type = 'prob')
        TMP_OBJECT[,19+i]<-as.numeric(iteration[,2])
        names(TMP_OBJECT)[19+i]<-paste0("p_TARGET_",i)
        a<-names(TMP_OBJECT)[19+i]
        ITERATION_FDR<-psmFDR(TMP_OBJECT, a, "DECOY",concat_decoy=0, prot_id="database")
        ITERATION_FDR_passed<-ITERATION_FDR[which(FIRST_ITERATION_FDR$psmFDR <=0.01),]
        df_temp<-ITERATION_FDR_passed[,c(1:19)]

        dataframes[[counter]]<-df_temp
        psms[counter]<-nrow(df_temp)
        peptides[counter]<-length(unique(paste(ITERATION_FDR_passed$Peptide)))
        proteins[counter]<-length(unique(paste(ITERATION_FDR_passed$Proteins)))

        cat("applying to full data", i,"\n")

        iteration_all<-predict(glmFit_2, dataPSMMat[,c(5,6,7,10,19)], type = 'prob')
        df[,19+i]<-as.numeric(iteration_all[,2])
        names(df)[19+i]<-paste0("p_TARGET_",i)
        a<-names(df)[19+i]
        ITERATION_FDR_all<-psmFDR(df, a, "DECOY",concat_decoy=0, prot_id="database")
        ITERATION_FDR_passed_all<-ITERATION_FDR_all[which(FIRST_ITERATION_FDR_all$psmFDR <=0.01),]
        df_temp_all<-ITERATION_FDR_passed_all[,c(1:19)]

        dataframes_all[[counter]]<-df_temp_all
        psms_all[counter]<-nrow(df_temp_all)
        peptides_all[counter]<-length(unique(paste(ITERATION_FDR_passed_all$Peptide)))
        proteins_all[counter]<-length(unique(paste(ITERATION_FDR_passed_all$Proteins)))

        only_targets<-df_temp_all[which(df_temp_all$Label=="TARGET"),]
        psms_all_clean[counter]<-nrow(only_targets)
        peptides_all_clean[counter]<-length(unique(paste(only_targets$Peptide)))
        proteins_all_clean[counter]<-length(unique(paste(only_targets$Proteins)))


        counter<-counter+1
    }
}

TRAINING_OBJECT$sum_na<-apply(TRAINING_OBJECT,1, FUN=function(x) sum(is.na(x)))

TRAINING_OBJECT<-TRAINING_OBJECT[which(TRAINING_OBJECT$sum_na==0),]
#nrow(606)

TRAINING_OBJECT$sum_na<-NULL
save(TRAINING_OBJECT, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/TRAINING_OBJECT.Rdata")

glmFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "glm", metric="ROC", trControl = ctrl)
iteration<-predict(glmFit, TRAINING_OBJECT[,c(3:6)], type = 'prob')  #primera

ALL_EXPERIMENT$sum_na<-apply(ALL_EXPERIMENT,1, FUN=function(x) sum(is.na(x)))
ALL_EXPERIMENT_COMPLETE_CASES<-ALL_EXPERIMENT[which(ALL_EXPERIMENT$sum_na==0),]

save(ALL_EXPERIMENT_COMPLETE_CASES, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/ALL_EXPERIMENT_COMPLETE_CASES.Rdata")

#### GLM
iteration<-predict(glmFit, ALL_EXPERIMENT_COMPLETE_CASES[,c(3:6)], type = 'prob')  #primera
ALL_EXPERIMENT_COMPLETE_CASES[,8]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT_COMPLETE_CASES)[8]<-"pred_TARGET"
glmroc<-roc(predictor=ALL_EXPERIMENT_COMPLETE_CASES$pred_TARGET, response=ALL_EXPERIMENT_COMPLETE_CASES$database, levels=rev(levels(ALL_EXPERIMENT_COMPLETE_CASES$database)))
plot(glmroc,main=paste0("glmROC AUC:", glmroc$auc))

#### rf
ctrl <- trainControl(method="oob", classProbs=TRUE, summaryFunction= twoClassSummary,savePredictions = TRUE,verboseIter=TRUE)


rfFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "rf",tuneLength=5, metric="ROC", trControl = ctrl,na.action = na.pass))
iteration<-predict(rfFit, ALL_EXPERIMENT_COMPLETE_CASES[,c(3:6)], type = 'prob')  #primera
ALL_EXPERIMENT_COMPLETE_CASES[,9]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT_COMPLETE_CASES)[9]<-"pred_TARGET_rf"

rfroc<-roc(predictor=ALL_EXPERIMENT_COMPLETE_CASES$pred_TARGET_rf, response=ALL_EXPERIMENT_COMPLETE_CASES$database, levels=rev(levels(ALL_EXPERIMENT_COMPLETE_CASES$database)))
plot(rfroc,main=paste0("rfROC AUC:",rfroc$auc))

#### BAYES
bayesFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "naive_bayes", metric="ROC", trControl = ctrl)
iteration<-predict(bayesFit, ALL_EXPERIMENT_COMPLETE_CASES[,c(3:6)], type = 'prob')  #primera
ALL_EXPERIMENT_COMPLETE_CASES[,10]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT_COMPLETE_CASES)[10]<-"pred_TARGET_bayes"

bayesroc<-roc(predictor=ALL_EXPERIMENT_COMPLETE_CASES$pred_TARGET_bayes, response=ALL_EXPERIMENT_COMPLETE_CASES$database, levels=rev(levels(ALL_EXPERIMENT_COMPLETE_CASES$database)))
plot(bayesroc,main=paste0("bayesROC AUC:",bayesroc$auc))

#### SVM
svmFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "svmRadial", metric="ROC", trControl = ctrl, preProcess = c("center", "scale"))
iteration<-predict(svmFit, ALL_EXPERIMENT_COMPLETE_CASES[,c(3:6)], type = 'prob')  #primera
ALL_EXPERIMENT_COMPLETE_CASES[,11]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT_COMPLETE_CASES)[11]<-"pred_TARGET_svm"

svmroc<-roc(predictor=ALL_EXPERIMENT_COMPLETE_CASES$pred_TARGET_svm, response=ALL_EXPERIMENT_COMPLETE_CASES$database, levels=rev(levels(ALL_EXPERIMENT_COMPLETE_CASES$database)))
plot(svmroc,main=paste0("svmROC AUC:",svmroc$auc))

### NNET

nnetGrid = expand.grid(size = seq(from = 1, to = 7, by = 1), decay = seq(from = 0.1, to = 0.5, by = 0.1))
nnetGrid = expand.grid(.size = c(5, 6, 7), .decay =c(0.5, 0.1))

NNETFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "nnet", metric="ROC", trControl = ctrl)
iteration<-predict(NNETFit, ALL_EXPERIMENT_COMPLETE_CASES[,c(3:6)], type = 'prob')  #primera
ALL_EXPERIMENT_COMPLETE_CASES[,12]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT_COMPLETE_CASES)[12]<-"pred_TARGET_NNET"

NNETroc<-roc(predictor=ALL_EXPERIMENT_COMPLETE_CASES$pred_TARGET_NNET, response=ALL_EXPERIMENT_COMPLETE_CASES$database, levels=rev(levels(ALL_EXPERIMENT_COMPLETE_CASES$database)))
plot(NNETroc,main=paste0("NNET_ROC AUC:",NNETroc$auc))

# "C5.0"

c5Fit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "C5.0", metric="ROC", trControl = ctrl, na.action = na.pass)
iteration<-predict(c5Fit, ALL_EXPERIMENT[,c(3:6)], type = 'prob', na.action = na.pass)  #primera
ALL_EXPERIMENT[,8]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT)[8]<-"pred_TARGET_c5"

c5roc<-roc(predictor=ALL_EXPERIMENT$pred_TARGET_c5, response=ALL_EXPERIMENT$database, levels=rev(levels(ALL_EXPERIMENT$database)))
plot(c5roc,main=paste0("c5_ROC AUC:",c5roc$auc))

### no va el predict
rfFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "rf",tuneLength=5, metric="ROC", trControl = ctrl,na.action = na.pass)
iteration<-predict(rfFit, ALL_EXPERIMENT[,c(3:6)], type = 'prob', na.action= na.pass)  #primera
ALL_EXPERIMENT_COMPLETE_CASES[,9]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT_COMPLETE_CASES)[9]<-"pred_TARGET_rf"

rfroc<-roc(predictor=ALL_EXPERIMENT_COMPLETE_CASES$pred_TARGET_rf, response=ALL_EXPERIMENT_COMPLETE_CASES$database, levels=rev(levels(ALL_EXPERIMENT_COMPLETE_CASES$database)))
plot(rfroc,main=paste0("rfROC AUC:",rfroc$auc))

#no va
glmFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "glm",tuneLength=5, metric="ROC", trControl = ctrl,na.action = na.pass)
iteration<-predict(glmFit, ALL_EXPERIMENT[,c(3:6)], type = 'prob', na.action= na.pass)  #primera
ALL_EXPERIMENT[,9]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT)[9]<-"pred_TARGET_glm"

# funciona pero no clasifica bien, será porque hay mucho NA?

bayesFit<-train(TRAINING_OBJECT[,c(3:6)], TRAINING_OBJECT[,7], method = "naive_bayes", metric="ROC", trControl = ctrl, na.action=na.pass)
iteration<-predict(bayesFit, ALL_EXPERIMENT[,c(3:6)], type = 'prob', na.action=na.pass)  #primera
ALL_EXPERIMENT[,9]<-as.numeric(iteration[,2])
names(ALL_EXPERIMENT)[9]<-"pred_TARGET_bayes"

bayesroc<-roc(predictor=ALL_EXPERIMENT$pred_TARGET_bayes, response=ALL_EXPERIMENT$database, levels=rev(levels(ALL_EXPERIMENT$database)))
plot(bayesroc,main=paste0("bayesROC AUC:",bayesroc$auc))

#### VOY A PROBAR CON 1 MOTOR, tandem)
> only_tandem<-TRAINING_OBJECT[which(!is.na(TRAINING_OBJECT$score_tandem)),]
> bayesFit<-train(only_tandem[,"score_tandem"], only_tandem[,"database"], method = "naive_bayes", metric="ROC", trControl = ctrl)
 iteration<-predict(bayesFit, ALL_EXPERIMENT[,c(3,6)], type = 'prob', laplace=0)  #primera


pdf("ROCs.pdf")
par(mfrow = c(3,3))
plot(glmroc,main=paste0("glmROC AUC:", glmroc$auc))
plot(rfroc,main=paste0("rfROC AUC:",rfroc$auc))
plot(bayesroc,main=paste0("bayesROC AUC:",bayesroc$auc))
plot(svmroc,main=paste0("svmROC AUC:",svmroc$auc))
plot(NNETroc,main=paste0("NNET_ROC AUC:",NNETroc$auc))
plot(c5roc,main=paste0("c5_ROC AUC:",c5roc$auc))

dev.off()
