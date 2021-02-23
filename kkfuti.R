



library(stringr)

data$Scanstring2<-gsub(" ","_",data$ScanString)
data$ScanNumber<-lapply(strsplit(paste(data$ScanString),"scan%3d"),"[",2)
data$ScanNumber<-str_extract(data$ScanNumber,"[0-9]+")
data$File<-lapply(strsplit(paste(data$Scanstring2),"\\."),"[",1)

data$Query<-paste(data$File,data$ScanNumber,sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq, sep="_")

data$ScanString<-data$Scanstring2
data$Scanstring2<-NULL
data$File<-NULL

save(data, file="results_Peptides_mascot_PXD003485_dataPSMMat.rda")
save(data, file="results_Peptides_mascot_PXD003485_dataPSMMat_filterAA.rda")
save(data, file="results_Peptides_mascot_PXD003485_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data, file="results_Peptides_mascot_PXD003485_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data, file="results_Peptides_mascot_PXD003485_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")

save(data, file="results_Peptides_mascot_PXD003560_dataPSMMat.rda")
save(data, file="results_Peptides_mascot_PXD003560_dataPSMMat_filterAA.rda")
save(data, file="results_Peptides_mascot_PXD003560_dataPSMMat_filterAA_PSMFDR_filter.rda")
save(data, file="results_Peptides_mascot_PXD003560_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter.rda")
save(data, file="results_Peptides_mascot_PXD003560_dataPSMMat_filterAA_psmFDR_Filter_pepFDR_Filter_protFDR_Filter.rda")


data$File<-lapply(strsplit(paste(data$ScanString),"\\."),"[",1)
data$File<-gsub("%2d","-",data$File)
data$File<-gsub("%","#",data$File)
data$Query<-paste(data$File,data$ScanNumber, sep="_")
data$File<-NULL


data$PSM<-gsub("%23","#", data$PSM)
data$Query<-gsub("%2d","-", data$Query)
data$Query<-gsub("%23","#",data$Query)

### 6271 mascot
data<-merge(data,dummydf, by="datfile", all.x=T)
data$mgf_file<-lapply(strsplit(paste(data$mgf_file),"\\."),"[",1)
data$Query<-paste(data$mgf_file,data$ScanString, sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq, sep="_")
data$mgf_file<-NULL


#omssa 6271
data$Query<-paste(lapply(strsplit(paste(data$datfile),"summary.csv"),"[",1),data$ScanString, sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq,sep="_")



tandemfiles<-c()
for(i in 1:length(experiments)){
    tandem_file<-paste0(tandem_folder, experiments[i],".txt")
    tandemfiles<-c(tandemfiles,tandem_file)
}

##xtandem results_Peptides_mascot_PXD003560_dataPSMMat

data$ScanString<-gsub(" ","_",data$ScanString)
data$ScanNumber<-lapply(strsplit(paste(data$ScanString),"\\."),"[",2)
data$Query<-paste(lapply(strsplit(paste(data$ScanString),"\\."),"[",1), data$ScanNumber, sep="_")
data$PSM<-paste(data$Query,data$PeptideSeq,sep="_")


aa<-list.files(getwd(), pattern="bestPSM.Rdata")


output_omssa<-data.frame()
for (i in 1:length(aa)){
    cat(aa[i],"\n")
    aa_loaded<-get(load(aa[i]))
    output_omssa<-rbind(output_omssa, aa_loaded)
}
