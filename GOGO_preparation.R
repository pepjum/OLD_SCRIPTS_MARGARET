

files<-list.files(pattern=".tsv")

df_rev<-data.frame()
output<-data.frame()

for(file in files){
    data_f<-read.table(file, header=T)
    n_proteins<-length(unique(paste(data_f$NX_code)))
    chr<-lapply(strsplit(paste(file),"_"),"[",3)
    datos<-data.frame("chr"=paste(chr), "proteins"=n_proteins)
    df_rev<-rbind(df_rev, datos)
    output<-rbind(output,data_f)
}

write.table(output, file="all_NX_GO.txt", col.names=T, quote=F, sep="\t")



##################################### GOGO preparation

all_NX<-read.table("~/data/pepe/21_neXtprot_2020_01_17/NEXTPROT_XML/xml/all_NX_GO.txt", header=T)


