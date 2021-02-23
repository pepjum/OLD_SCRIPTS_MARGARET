args=(commandArgs(TRUE))

library(stringr)

ALL_EXPERIMENT_file<-args[1]
#

#batch_number<-str_extract(tandem_file, "[0-9]+")


ALL_EXPERIMENT_file_batch<-get(load(ALL_EXPERIMENT_file))

output_databases<-data.frame()
output_PeptideMASS<-data.frame()
dummys_peptides<-c()
dummys_databases<-c()
for(row in 1:nrow(ALL_EXPERIMENT_file_batch)){
    tmp_row<-ALL_EXPERIMENT_file_batch[row,]
    for(i in c(2,5,8)){
        print(tmp_row[,i])
        if(is.na(tmp_row[,i])){
        next
        }else{
            if (tmp_row$dummy %in% dummys_peptides){next}else{
            print(length(dummys_peptides))
            PeptideMass<-tmp_row[,i]
            print(PeptideMass)
            dummys_peptides<-c(dummys_peptides,tmp_row$dummy)
            out<-data.frame("PeptideMass"=PeptideMass,"dummy"=tmp_row$dummy)
            cat(paste(out),"\n")
            output_PeptideMASS<-rbind(output_PeptideMASS,out)
        }
        }
    }
    for(i in c(4,7,10)){
        if(tmp_row[,i]=="NA"){
        next
        }else{
            if (tmp_row$dummy %in% dummys_databases){next}else{

            database<-tmp_row[,i]
            dummys_databases<-c(dummys_databases,tmp_row$dummy)
            out_db<-data.frame("database"=database,"dummy"=tmp_row$dummy)
            #cat(paste(out),"\n")
            output_databases<-rbind(output_databases,out_db)
        }
        }
    }
}

out_out<-merge(output_databases, output_PeptideMASS, by="dummy")
out_out$dummy<-as.numeric(out_out$dummy)
if(nrow(out_out) != nrow(ALL_EXPERIMENT_file_batch)){
    output_file<-paste0(lapply(strsplit(paste(ALL_EXPERIMENT_file),"\\."),"[",1),"_log.txt")
    a<-paste0("something is wrong in ",ALL_EXPERIMENT_file)
    write.table(a, output_file)
}

output_file<-paste0(lapply(strsplit(paste(ALL_EXPERIMENT_file),"\\."),"[",1),"_PeptideMass_AND_DB.Rdata")
save(out_out, file=output_file)
