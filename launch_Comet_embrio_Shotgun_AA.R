

library(doParallel)

cl <- makeCluster(5)

Experiments<-c("PXD004003","PXD004246","PXD004288","PXD004494","PXD004938")

files<-paste0("Shotgun_comet_",Experiments,".txt")
files<-paste("/home/margaret/data/pepe/13_EMBRIO_NOV18/Comet_files",files,sep="/")

for(j in 1:length(files)){
        y<- paste0('Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/ShotgunAnalysisAutomated_margaret.R ', files[j]," ", Experiments[j],' 0'," D_sp"," 3"," 2"," 9"," 50")
        cat(y,"\n")
        system(y)

}
stopCluster(cl)



#### mascot


Directories<-list.dirs()
Directories<-grep("*-D", Directories, value = TRUE)
dirs<-Directories[-c(43:51)]
dirs<-lapply(strsplit(paste(dirs),"./"),"[", 2)
dirs<-lapply(strsplit(paste(dirs),"-D"),"[", 1)

files<-paste("/home/margaret/data/pepe/13_EMBRIO_NOV18/Dat_Files",dirs,sep="/")
Experiments<-basename(files)

output<-c()

for(j in 1:length(files)){
        y<- paste0('Rscript /home/margaret/data/pepe/scripts/ANALISIS_PROTEOMICA/ShotgunAnalysisAutomated_margaret.R ', "/home/margaret/data/pepe/13_EMBRIO_NOV18/Dat_Files/all_files.txt"," ", Experiments[j],' 0'," DECOY"," 3"," 1"," 9"," 50", " &")
        #cat(y,"\n")
        #system(y)
        output<-c(output,y)

}

write.table(output, file="/home/margaret/data/pepe/13_EMBRIO_NOV18/Dat_Files/launch_SAA.txt", col.names=F,row.names=F,sep="\t",quote=F)
