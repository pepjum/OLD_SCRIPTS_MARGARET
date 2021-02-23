
mascot_best_psm<-get(load("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/All_mascot_best_PSM_NCI60_DP.Rdata"))
comet_best_psm<-get(load("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/All_comet_best_PSM_NCI60_DP.Rdata"))
tandem_best_psm<-get(load("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/All_tandem_best_PSM_NCI60_DP.Rdata"))
omssa_best_psm<-get(load("/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/All_omssa_best_PSM_NCI60_DP.Rdata"))




unique_mascot_PSM<-unique(paste(mascot_best_psm$PSM))
unique_comet_PSM<-unique(paste(comet_best_psm$PSM))
unique_tandem_PSM<-unique(paste(tandem_best_psm$PSM))
unique_omssa_PSM<-unique(paste(omssa_best_psm$PSM))


list_querys<-c(unique_mascot_PSM,unique_comet_PSM, unique_tandem_PSM, unique_omssa_PSM)
list_querys_unique<-unique(list_querys)  #348152

mascot_tmp<-mascot_best_psm[,c(21,6,8,10,16)]
names(mascot_tmp)<-c("PSM","PeptideMass_mascot","PeptideSeq_mascot","score_mascot","database_mascot")
mascot_tmp<-unique(mascot_tmp)
mascot_tmp$PeptideMass_mascot<-as.numeric(paste(mascot_tmp$PeptideMass_mascot))
mascot_tmp$score_mascot<-as.numeric(paste(mascot_tmp$score_mascot))
mascot_tmp$PeptideSeq_mascot<-as.character(mascot_tmp$PeptideSeq_mascot)
mascot_tmp$database_mascot<-as.character(mascot_tmp$database_mascot)
#mascot_tmp$psmFDR_mascot<-psmFDR(mascot_tmp, "score_mascot", "DECOY", 0, "database_mascot")

comet_tmp<-comet_best_psm[,c(26,8,10,12,18)]
names(comet_tmp)<-c("PSM","PeptideMass_comet","PeptideSeq_comet","score_comet","database_comet")
comet_tmp<-unique(comet_tmp)
comet_tmp$PeptideMass_comet<-as.numeric(paste(comet_tmp$PeptideMass_comet))
comet_tmp$score_comet<-as.numeric(paste(comet_tmp$score_comet))
comet_tmp$PeptideSeq_comet<-as.character(comet_tmp$PeptideSeq_comet)
comet_tmp$database_comet<-as.character(comet_tmp$database_comet)
#comet_tmp$psmFDR_comet<-psmFDR(comet_tmp, "score_comet", "DECOY", 0, "database_comet")

tandem_tmp<-tandem_best_psm[,c(23,7,9,11,18)]
names(tandem_tmp)<-c("PSM","PeptideMass_tandem","PeptideSeq_tandem","score_tandem","database_tandem")
tandem_tmp<-unique(tandem_tmp)
#tandem_tmp$psmFDR_tandem<-psmFDR(comet_tmp, "score_tandem", "DECOY", 0, "database_tandem")
tandem_tmp$PeptideMass_tandem<-as.numeric(paste(tandem_tmp$PeptideMass_tandem))
tandem_tmp$score_tandem<-as.numeric(paste(tandem_tmp$score_tandem))
tandem_tmp$PeptideSeq_tandem<-as.character(tandem_tmp$PeptideSeq_tandem)
tandem_tmp$database_tandem<-as.character(tandem_tmp$database_tandem)


omssa_tmp<-omssa_best_psm[,c(22,7,9,11,17)]
names(omssa_tmp)<-c("PSM","PeptideMass_omssa","PeptideSeq_omssa","score_omssa","database_omssa")
omssa_tmp<-unique(omssa_tmp)

omssa_tmp$PeptideMass_omssa<-as.numeric(paste(omssa_tmp$PeptideMass_omssa))
omssa_tmp$score_omssa<-as.numeric(paste(omssa_tmp$score_omssa))
omssa_tmp$PeptideSeq_omssa<-as.character(omssa_tmp$PeptideSeq_omssa)
omssa_tmp$database_omssa<-as.character(omssa_tmp$database_omssa)




Query_df<-data.frame("PSM"=list_querys_unique)

Query_df<-merge(Query_df, mascot_tmp, by.x="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, comet_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, tandem_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)
Query_df<-merge(Query_df, omssa_tmp,by="PSM", all.x=T)
Query_df<-unique(Query_df)

save(Query_df, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/MATRIX_PSM_initial.Rdata")


Query_df$PeptideSeq<-lapply(strsplit(paste(Query_df$PSM),"_"),"[",8)

Query_df$PeptideSeq_comet<-NULL
Query_df$PeptideSeq_tandem<-NULL
Query_df$PeptideSeq_mascot<-NULL
Query_df$PeptideSeq_omssa<-NULL



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

Query_df$database<-"NA"
for(row in 1:nrow(Query_df)){
    cat(row,"\n")
    b<-c(Query_df$database_mascot[row], Query_df$database_comet[row], Query_df$database_omssa[row], Query_df$database_tandem[row])
    c<-na.omit(b)[1]
    Query_df$database[row]<-c
}

Query_df$database_mascot<-NULL
Query_df$database_comet<-NULL
Query_df$database_omssa<-NULL
Query_df$database_tandem<-NULL


save(Query_df, file="/home/margaret/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/MATRIZ_BEST_PSM_NCI60_DP.Rdata")
