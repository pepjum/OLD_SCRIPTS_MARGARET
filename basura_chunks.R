results_mascot$Query<-gsub("%2d","-",results_mascot$Query)
results_mascot$Query<-gsub("-","_",results_mascot$Query)
results_mascot$Query<-gsub("#","_",results_mascot$Query)


results_mascot$Query<-toupper(results_mascot$Query)

results_comet$Query<-paste(results_comet$Fraction,results_comet$ScanNumber,sep="_")
results_comet$Query<-toupper(results_comet$Query)
results_comet$Query<-gsub("-","_",results_comet$Query)
results_comet$Query<-gsub("#","_",results_comet$Query)



results_omssa$Query<-toupper(results_omssa$Query)
results_omssa$Query<-gsub("-","_",results_omssa$Query)
results_omssa$Query<-gsub("#","_",results_omssa$Query)


results_tandem$Query<-toupper(results_tandem$Query)
results_tandem$Query<-gsub("-","_",results_tandem$Query)
results_tandem$Query<-gsub("#","_",results_tandem$Query)



results_omssa$tmp<-lapply(strsplit(paste(results_omssa$datfile),"_summary"),"[",1)
results_omssa$PeptideSeq<-toupper(results_omssa$PeptideSeq)
results_omssa$Query<-paste(results_omssa$tmp,results_omssa$ScanString,sep="_")
results_omssa$PSM<-paste(results_omssa$Query,results_omssa$PeptideSeq,sep="_")
results_omssa$tmp<-NULL

all<-merge(data,trazability, by="datfile", all.x=T)
all$tmp<-lapply(strsplit(paste(all$mgf_file)),".mg"),"[",1)
all$Query<-paste(all$tmp, all$ScanString, sep="_")
all$PSM<-paste(all$Query, all$PeptideSeq, sep="_")
all$tmp<-NULL

save(all, file="/home/margaret/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD006271_dataPSMMat.rda")

results_tandem$ScanString<-lapply(strsplit(paste(results_tandem$ScanString),"_"),"[",1)
results_tandem$Query<-paste(lapply(strsplit(paste(results_tandem$datfile),"\\."),"[",1),results_tandem$ScanString,sep="_")
results_tandem$PSM<-paste(results_tandem$Query,results_tandem$PeptideSeq,sep="_")

save(results_tandem,file="/home/margaret/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD006271_dataPSMMat.rda")
save(results_tandem,file="/home/margaret/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD006271_dataPSMMat_filterAA.rda")
save(results_tandem,file="/home/margaret/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD006271_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(results_tandem,file="/home/margaret/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD006271_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(results_tandem,file="/home/margaret/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD006271_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")

###PXD006271
data_merged<-merge(data,trazability,by.x="datfile",by.y="datFile", all.x=T)
data<-data_merged
data$temp<-lapply(strsplit(paste(data$MGFFile),"\\."),"[",1)
data$PSM<-paste(data$temp, data$ScanString, data$PeptideSeq, sep="_")
data$Query<-paste(data$temp, data$ScanString, sep="_")
data$temp<-NULL

save(data, file="results_Peptides_mascot_PXD006271_dataPSMMat.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD006271_dataPSMMat_filterAA.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD006271_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD006271_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda"")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD006271_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")

data$PSM<-str_replace_all(data$PSM,"%23","#")


###PXD009893_Q_Exactive Mascot

library("readxl")

excel<-read_excel("13_EMBRIO_NO18_t2.xls")
trazability<-excel[,c("datFile","MGFFile")]
trazability<-as.data.frame(trazability)

data$ScanNumber<-lapply(strsplit(paste(data$ScanString)," "),"[", 7)
data$ScanNumber<-str_extract(data$ScanNumber,"[0-9]+")
data_merged<-merge(data,trazability,by.x="datfile",by.y="datFile", all.x=T)
data<-data_merged
data$temp<-lapply(strsplit(paste(data$MGFFile),"\\."),"[",1)
data$PSM<-paste(data$temp, data$ScanNumber, data$PeptideSeq, sep="_")
data$Query<-paste(data$temp, data$ScanNumber, sep="_")
data$temp<-NULL

save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_Q_Exactive_dataPSMMat.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_Q_Exactive_dataPSMMat_filterAA.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_Q_Exactive_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_Q_Exactive_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_Q_Exactive_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter-protFDR_Filter.rda")

data$ScanNumber<-lapply(strsplit(paste(data$ScanString)," "),"[", 1)
data_merged<-merge(data,trazability,by.x="datfile",by.y="datFile", all.x=T)
data<-data_merged
data$temp<-lapply(strsplit(paste(data$MGFFile),"\\."),"[",1)
data$PSM<-paste(data$temp, data$ScanNumber, data$PeptideSeq, sep="_")
data$Query<-paste(data$temp, data$ScanNumber, sep="_")
data$temp<-NULL

save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_T_TOF_dataPSMMat.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_T_TOF_dataPSMMat_filterAA.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_T_TOF_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_T_TOF_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD009893_T_TOF_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")


data$Query<-paste(lapply(strsplit(paste(data$datfile),"\\."),"[",1),lapply(strsplit(paste(data$ScanString),"_"),"[",1),sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq,sep="_")

save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_T_TOF_dataPSMMat.rda")
save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_T_TOF_dataPSMMat_filterAA.rda")
save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_T_TOF_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_T_TOF_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_T_TOF_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")


data$tmp<-lapply(strsplit(paste(data$ScanString),";_"),"[",3)
data$ScanNumber<-lapply(strsplit(paste(data$tmp),"_"),"[",2)
data$Query<-paste(lapply(strsplit(paste(data$datfile),"\\."),"[",1),data$ScanNumber,sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq,sep="_")
data$tmp<-NULL

save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_Q_Exactive_dataPSMMat.rda")
save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_Q_Exactive_dataPSMMat_filterAA.rda")
save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_Q_Exactive_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_Q_Exactive_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data,file="~/data/pepe/13_EMBRIO_NOV18/Tandem_Files/results_Peptides_tandem_PXD009893_Q_Exactive_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")


data$ScanString<-gsub(" ","_",data$ScanString)
data$ScanNumber<-lapply(strsplit(paste(data$ScanString),"\\."),"[",2)
data$Query<-paste(lapply(strsplit(paste(data$ScanString),"\\."),"[",1),lapply(strsplit(paste(data$ScanString),"\\."),"[",2),sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq,sep="_")

save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003560_dataPSMMat.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003560_dataPSMMat_filterAA.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003560_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003560_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003560_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")

save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003485_dataPSMMat.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003485_dataPSMMat_filterAA.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003485_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003485_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data, file="~/data/pepe/13_EMBRIO_NOV18/Dat_Files/results_Peptides_mascot_PXD003485_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")


data$ScanNumber<-lapply(strsplit(paste(data$scanstring),"\\."),"[",2)
data$Query<-paste(data$Fraction,data$ScanNumber,sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq,sep="_")

data$ScanNumber<-lapply(strsplit(paste(data$ScanString),"\\."),"[",2)
data$Query<-paste(lapply(strsplit(paste(data$ScanString),"\\."),"[",1), data$ScanNumber, sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq,sep="_")



save(data,file="results_Peptides_tandem_PXD005123_dataPSMMat.rda")
save(data,file="results_Peptides_tandem_PXD005123_dataPSMMat_filterAA.rda")
save(data,file="results_Peptides_tandem_PXD005123_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data,file="results_Peptides_tandem_PXD005123_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data,file="results_Peptides_tandem_PXD005123_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")



files<-files<-list.files(pattern="_bestPSM.Rdata")

output<-data.frame()
for(i in 1:length(files)){
    cat(files[i],"\n")
    data<-get(load(files[i]))
    output<-rbind(output,data)
}

save(output, file="All_comet_best_PSM_NCI60_DP.Rdata")
save(output, file="All_omssa_best_PSM_NCI60_DP.Rdata")
save(output, file="All_tandem_best_PSM_NCI60_DP.Rdata")
save(output, file="All_mascot_best_PSM_NCI60_DP.Rdata")

save(output, file="All_comet_NCI60_DP.Rdata")
save(output, file="All_omssa_NCI60_DP.Rdata")
save(output, file="All_tandem_NCI60_DP.Rdata")
save(output, file="All_mascot_NCI60_DP.Rdata")
