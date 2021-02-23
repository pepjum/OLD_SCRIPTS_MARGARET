
mascot_raw<-get(load("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_BEST_PSM_MASCOT.Rdata"))
comet_raw<-get(load("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_BEST_PSM_COMET.Rdata"))
omssa_raw<-get(load("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_BEST_PSM_OMSSA.Rdata"))
tandem_raw<-get(load("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_BEST_PSM_TANDEM.Rdata"))

source("/home/margaret/data/pepe/scripts/functions_Classifier.R")
library(stringr)
library(plyr)
library(dplyr)


unique_mascot_PSM<-unique(paste(mascot_raw$PSM))
unique_comet_PSM<-unique(paste(comet_raw$PSM))
unique_omssa_PSM<-unique(paste(omssa_raw$PSM))
unique_tandem_PSM<-unique(paste(tandem_raw$PSM))

list_querys<-c(unique_mascot_PSM,unique_comet_PSM,unique_omssa_PSM, unique_tandem_PSM)
list_querys_unique<-unique(list_querys)   #1285431


mascot_tmp<-mascot_raw[,c(19,4,8,14)]
names(mascot_tmp)<-c("PSM","PeptideMass_mascot","score_mascot","database_mascot")
mascot_tmp<-unique(mascot_tmp)
#mascot_tmp$psmFDR_mascot<-psmFDR(mascot_tmp, "score_mascot", "DECOY", 0, "database_mascot")

comet_tmp<-comet_raw[,c(24,6,10,16)]
names(comet_tmp)<-c("PSM","PeptideMass_comet","score_comet","database_comet")
comet_tmp<-unique(comet_tmp)
#comet_tmp$psmFDR_comet<-psmFDR(comet_tmp, "score_comet", "DECOY", 0, "database_comet")

tandem_tmp<-tandem_raw[,c(21,5,15,16)]
names(tandem_tmp)<-c("PSM","PeptideMass_tandem","score_tandem","database_tandem")
tandem_tmp<-unique(tandem_tmp)

omssa_tmp<-omssa_raw[,c(20,5,9,15)]
names(omssa_tmp)<-c("PSM","PeptideMass_omssa","score_omssa","database_omssa")
omssa_tmp<-unique(omssa_tmp)

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

Query_df$database_mascot<-as.character(paste(Query_df$database_mascot))
Query_df$database_comet<-as.character(paste(Query_df$database_comet))
Query_df$database_tandem<-as.character(paste(Query_df$database_tandem))
Query_df$database_omssa<-as.character(paste(Query_df$database_omssa))

for(number in c(2,3,5,6,8,9,11,12)){
    cat(number,"\n")
    Query_df[,number]<-as.numeric(paste(Query_df[,number]))
}

Query_df$PeptideMass<-"NA"
for(row in 1:nrow(Query_df)){
    cat(row,"\n")
    b<-c(Query_df$PeptideMass_mascot[row], Query_df$PeptideMass_comet[row], Query_df$PeptideMass_tandem[row],Query_df$PeptideMass_omssa[row])
    c<-na.omit(b)[1]
    Query_df$PeptideMass[row]<-c
}


Query_df$PeptideMass_mascot<-NULL
Query_df$PeptideMass_comet<-NULL
Query_df$PeptideMass_omssa<-NULL
Query_df$PeptideMass_tandem<-NULL

save(Query_df, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/MATRIX_PSM_initial.Rdata")


cat("preparing training set....\n")

ALL_DECOY<-rbind(Query_df[which(Query_df$database_mascot=="D"),],Query_df[which(Query_df$database_comet=="D"),],Query_df[which(Query_df$database_tandem=="D"),], Query_df[which(Query_df$database_omssa=="D"),])
ALL_DECOY<-unique(ALL_DECOY)
ALL_DECOY<-ALL_DECOY[!(ALL_DECOY$database_mascot=="T"),]
ALL_DECOY<-ALL_DECOY[!(ALL_DECOY$database_comet=="T"),]
ALL_DECOY<-ALL_DECOY[!(ALL_DECOY$database_omssa=="T"),]
ALL_DECOY<-ALL_DECOY[!(ALL_DECOY$database_tandem=="T"),]

save(ALL_DECOY, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_DECOY.Rdata")

ALL_DECOY <- mutate_all(ALL_DECOY, funs(replace(., .=='NA', NA)))

ALL_DECOY$database<-"NA"
for(row in 1:nrow(ALL_DECOY)){
    cat(row,"\n")
    b<-c(ALL_DECOY$database_mascot[row], ALL_DECOY$database_comet[row], ALL_DECOY$database_tandem[row],ALL_DECOY$database_omssa[row])
    c<-na.omit(b)[1]
    ALL_DECOY$database[row]<-c
}


ALL_DECOY$database_mascot<-NULL
ALL_DECOY$database_comet<-NULL
ALL_DECOY$database_omssa<-NULL
ALL_DECOY$database_tandem<-NULL


save(ALL_DECOY, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_DECOY_PREPARED.Rdata")



DECOYS_for_TRAINING<-sample_n(ALL_DECOY, nrow(ALL_DECOY)/2)
# DECOYS_for_TRAINING <- mutate_all(DECOYS_for_TRAINING, funs(replace(., .=='NA', NA)))
#
# DECOYS_for_TRAINING$database<-"NA"
# for(row in 1:nrow(DECOYS_for_TRAINING)){
#     cat(row,"\n")
#     b<-c(DECOYS_for_TRAINING$database_mascot[row], DECOYS_for_TRAINING$database_comet[row], DECOYS_for_TRAINING$database_tandem[row],DECOYS_for_TRAINING$database_omssa[row])
#     c<-na.omit(b)[1]
#     DECOYS_for_TRAINING$database[row]<-c
# }
#
#
# DECOYS_for_TRAINING$database_mascot<-NULL
# DECOYS_for_TRAINING$database_comet<-NULL
# DECOYS_for_TRAINING$database_omssa<-NULL
# DECOYS_for_TRAINING$database_tandem<-NULL

save(DECOYS_for_TRAINING, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/DECOYS_FOR_TRAINING_SET_FOR_FIRST_ITERATION_PREPARED.Rdata")




ALL_TARGET<-rbind(Query_df[which(Query_df$database_mascot=="T"),], Query_df[which(Query_df$database_comet=="T"),], Query_df[which(Query_df$database_tandem=="T"),], Query_df[which(Query_df$database_omssa=="T"),] )
ALL_TARGET<-unique(ALL_TARGET)
ALL_TARGET<-ALL_TARGET[!(ALL_TARGET$database_mascot=="D"),]
ALL_TARGET<-ALL_TARGET[!(ALL_TARGET$database_comet=="D"),]
ALL_TARGET<-ALL_TARGET[!(ALL_TARGET$database_omssa=="D"),]
ALL_TARGET<-ALL_TARGET[!(ALL_TARGET$database_tandem=="D"),]



save(ALL_TARGET, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_TARGET.Rdata")

ALL_TARGET <- mutate_all(ALL_TARGET, funs(replace(., .=='NA', NA)))

ALL_TARGET$database<-"NA"
for(row in 1:nrow(ALL_TARGET)){
    cat(row,"\n")
    b<-c(ALL_TARGET$database_mascot[row], ALL_TARGET$database_comet[row], ALL_TARGET$database_tandem[row],ALL_TARGET$database_omssa[row])
    c<-na.omit(b)[1]
    ALL_TARGET$database[row]<-c
}


ALL_TARGET$database_mascot<-NULL
ALL_TARGET$database_comet<-NULL
ALL_TARGET$database_omssa<-NULL
ALL_TARGET$database_tandem<-NULL

save(ALL_TARGET, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_TARGET_PREPARED.Rdata")

ALL_EXPERIMENT_NO_AMBIGUOUS_PSM<-rbind(ALL_TARGET,ALL_DECOY)

save(ALL_EXPERIMENT_NO_AMBIGUOUS_PSM, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_EXPERIMENT_WITHOUT_AMBIGUOUS_TyD_PREPARED.Rdata")



TMP_OBJECT<-rbind(ALL_TARGET, DECOYS_for_TRAINING)    #ESTE objecto es sobre el que hay que hacer primer FDR
TMP_OBJECT<-unique(TMP_OBJECT)    #nrow(TMP_OBJECT)

save(TMP_OBJECT, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/Rerank_object_prepared.Rdata")
### subset de este objeto para cada motor para poder calcular el FDR

TMP_OBJECT_mascot_selected<-TMP_OBJECT[which(!is.na(TMP_OBJECT$score_mascot)),]
mascot_DF_with_FDR<-psmFDR(TMP_OBJECT_mascot_selected, "score_mascot", "D", concat_decoy=0, prot_id="database")
mascot_TMP_FDR<-mascot_DF_with_FDR[,c("PSM","psmFDR")]
names(mascot_TMP_FDR)[2]<-"psmFDR_mascot"

TMP_OBJECT<-left_join(TMP_OBJECT, mascot_TMP_FDR, by = "PSM")

# TMP_OBJECT<-merge(TMP_OBJECT, mascot_TMP_FDR, by="PSM", all.x=TRUE)
TMP_OBJECT<-unique(TMP_OBJECT)

TMP_OBJECT_comet_selected<-TMP_OBJECT[which(!is.na(TMP_OBJECT$score_comet)),]
comet_DF_with_FDR<-psmFDR(TMP_OBJECT_comet_selected, "score_comet", "D", concat_decoy=0, prot_id="database")
comet_TMP_FDR<-comet_DF_with_FDR[,c("PSM","psmFDR")]
names(comet_TMP_FDR)[2]<-"psmFDR_comet"

TMP_OBJECT<-merge(TMP_OBJECT, comet_TMP_FDR, by="PSM", all.x=TRUE)

TMP_OBJECT<-unique(TMP_OBJECT)

TMP_OBJECT_tandem_selected<-TMP_OBJECT[which(!is.na(TMP_OBJECT$score_tandem)),]
tandem_DF_with_FDR<-psmFDR(TMP_OBJECT_tandem_selected, "score_tandem", "D", concat_decoy=0, prot_id="database")
tandem_TMP_FDR<-tandem_DF_with_FDR[,c("PSM","psmFDR")]
names(tandem_TMP_FDR)[2]<-"psmFDR_tandem"

TMP_OBJECT<-merge(TMP_OBJECT, tandem_TMP_FDR, by="PSM", all.x=TRUE)
TMP_OBJECT<-unique(TMP_OBJECT)

TMP_OBJECT_omssa_selected<-TMP_OBJECT[which(!is.na(TMP_OBJECT$score_omssa)),]
omssa_DF_with_FDR<-psmFDR(TMP_OBJECT_omssa_selected, "score_omssa", "D", concat_decoy=0, prot_id="database")
omssa_TMP_FDR<-omssa_DF_with_FDR[,c("PSM","psmFDR")]
names(omssa_TMP_FDR)[2]<-"psmFDR_omssa"

TMP_OBJECT<-merge(TMP_OBJECT, omssa_TMP_FDR, by="PSM", all.x=TRUE)
TMP_OBJECT<-unique(TMP_OBJECT)

save(TMP_OBJECT, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/TRAINING_SET_FOR_FIRST_ITERATION_With_FDR.Rdata")

rm(mascot_TMP_FDR,mascot_DF_with_FDR,comet_DF_with_FDR,comet_TMP_FDR,tandem_DF_with_FDR,tandem_TMP_FDR, omssa_TMP_FDR)

DF_FDR_PASSED_mascot<-TMP_OBJECT[which(na.omit(paste(TMP_OBJECT$psmFDR_mascot))<=0.01),]
DF_FDR_PASSED_comet<-TMP_OBJECT[which(na.omit(paste(TMP_OBJECT$psmFDR_comet))<=0.01),]
DF_FDR_PASSED_tandem<-TMP_OBJECT[which(na.omit(paste(TMP_OBJECT$psmFDR_tandem))<=0.01),]
DF_FDR_PASSED_omssa<-TMP_OBJECT[which(na.omit(paste(TMP_OBJECT$psmFDR_omssa))<=0.01),]


DF_FDR_PASSED<-unique(rbind(DF_FDR_PASSED_comet,DF_FDR_PASSED_mascot,DF_FDR_PASSED_tandem, DF_FDR_PASSED_omssa))

save(DF_FDR_PASSED, file="/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/TRAINING_SET_FOR_FIRST_ITERATION_PSMFDR_filtered.Rdata")   # con este objeto + mitad decoys se entrena en cada iteracion, rerankeandolo en cada iteracion


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
TRAINING_OBJECT$psmFDR_omssa<-NULL

TRAINING_OBJECT<-rbind(TRAINING_OBJECT,DECOYS_for_TRAINING)
save(TRAINING_OBJECT,file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/TRAINING_SET_FOR_FIRST_ITERATION.Rdata")

#### QUitar las redundancias de las columnas databases y dejar solo una con esa informacion
