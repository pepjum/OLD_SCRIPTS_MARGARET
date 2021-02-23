### uPes automatizer
library(optparse)
library(icesTAF)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Query object", metavar="character"),
  #make_option(c("-u", "--uPE"), type="character", default=NULL,
  #            help="uPE object", metavar="character"),
  #make_option(c("-p", "--PE1"), type="character", default=NULL,
  #            help="PE1 without uPE", metavar="character"),
  make_option(c("-d", "--database"), type="character", default=NULL,
              help="database TGCA, GTEX.... whatever", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

load(opt$file)
#upe1<-get(load(opt$uPE))
#pe1_no_upe<-get(load(opt$PE1))
# Primero voy a generar TCGA_cor, TCGA_corP y TCGA_corFDR, que los necesito para el analisis. Son lineas de codigo de ELI
#ensg_upe1<-paste(upe1$ENSG)
#ensg_pe1_sin_upe1<-paste(pe1_no_upe$ENSG)

output_dir<-paste0(dirname(opt$file), "/",opt$database,"/")
dir.create(output_dir, showWarnings = FALSE)

if(opt$database =="TCGA"){
    ensg_upe1<-paste(ensg_upe1_tcga)
    ensg_pe1_sin_upe1<-paste(ensg_pe1_sin_upe1_tcga)


    TCGA_sel_upe1<-intersect(ensg_upe1, rownames(TCGA_data))
    TCGA_sel_PE1<-intersect(ensg_pe1_sin_upe1, rownames(TCGA_data))
    TCGA_cor <- apply(TCGA_data[TCGA_sel_upe1,], 1, FUN=function(x) apply(TCGA_data[TCGA_sel_PE1,], 1, FUN=function(y) cor(x,y)))
    # write.table(rownames(TCGA_cor), file=paste0(output_dir,"TCGA_PE_names.txt"), row.names=F, col.names=F, quote=F, sep="\n")
    # write.table(colnames(TCGA_cor), file=paste0(output_dir,"TCGA_UPE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")
    save(TCGA_cor, file=paste0(output_dir,"TCGA_cor",".Rdata"))
    TCGA_corP <- apply(TCGA_data[TCGA_sel_upe1,], 1, FUN=function(x) apply(TCGA_data[TCGA_sel_PE1,], 1, FUN=function(y) cor.test(x,y)$p.value))
    save(TCGA_corP, file=paste0(output_dir,"TCGA_corP",".Rdata"))
    TCGA_corFDR <- apply(TCGA_corP, 2, FUN=function(x) p.adjust(x,method="fdr"))
    save(TCGA_corFDR, file=paste0(output_dir,"TCGA_corFDR",".Rdata"))
    TCGA_c <- (apply(TCGA_corFDR, 2, FUN=function(x) (x<0.05)*1)) * (apply(TCGA_cor, 2, FUN=function(x) (abs(x)>0.3)*1)) # meanGS=1081.105, medianGS=851  #eli
    save(TCGA_c, file=paste0(output_dir,"TCGA_c",".Rdata"))
    write.table(rownames(TCGA_c), file=paste0(output_dir,"TCGA_PE_names.txt"), row.names=F, col.names=F, quote=F, sep="\n")
    write.table(colnames(TCGA_c), file=paste0(output_dir,"TCGA_UPE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")

                                                                                                                        #summary(colSums(TCGA_c))     # primera
                                                                                                                        #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
                                                                                                                        #    0     289     752    1007    1525    4166
                                                                                                                        #nuevo summary(colSums(TCGA_c))
                                                                                                                        #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
                                                                                                                        #       3     586    1128    1287    1874    4308 

}else if(opt$database =="GTEX"){
    ensg_upe1<-paste(ensg_upe1_gtex)
    ensg_pe1_sin_upe1<-paste(ensg_pe1_sin_upe1_gtex)

    GTEX_sel_upe1<-intersect(ensg_upe1, rownames(GTEX_data))
    GTEX_sel_PE1<-intersect(ensg_pe1_sin_upe1, rownames(GTEX_data))

    GTEX_cor <- apply(GTEX_data[GTEX_sel_upe1,], 1, FUN=function(x) apply(GTEX_data[GTEX_sel_PE1,], 1, FUN=function(y) cor(x,y)))
    save(GTEX_cor, file=paste0(output_dir,"GTEX_cor",".Rdata"))
    write.table(rownames(GTEX_cor), file=paste0(output_dir,"GTEX_PE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")
    write.table(colnames(GTEX_cor), file=paste0(output_dir,"GTEX_UPE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")

    GTEX_corP <- apply(GTEX_data[GTEX_sel_upe1,], 1, FUN=function(x) apply(GTEX_data[GTEX_sel_PE1,], 1, FUN=function(y) cor.test(x,y)$p.value))
    save(GTEX_corP, file=paste0(output_dir,"GTEX_corP",".Rdata"))
    GTEX_corFDR <- apply(GTEX_corP, 2, FUN=function(x) p.adjust(x,method="fdr"))
    save(GTEX_corFDR, file=paste0(output_dir,"GTEX_corFDR",".Rdata"))
    GTEX_c <- (apply(GTEX_corFDR, 2, FUN=function(x) (x<0.05)*1)) * (apply(GTEX_cor, 2, FUN=function(x) (abs(x)>0.4)*1)) # meanGS=1151, medianGS=1336.135, 0_GS=3. #eli
                                                                                                                         #meanGS=2905, medianGS=2905   0_GS= 3
                                                                                                                         #subo a 0.4 summary(colSums(GTEX_c))
                                                                                                                        # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   0_Gs=15
                                                                                                                        #0.0   424.8  1027.0  1264.1  1929.2  4360.0

                                                                                                                        #nuevo 0.3
                                                                                                                        #summary(colSums(dat))    0_GS=0
                                                                                                                        #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
                                                                                                                        #166    2130    3013    3064    4016    6275 
                                                                                                                        #0.4
#                                                                                                                         summary(colSums(GTEX_c))
                                                                                                                            #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
                                                                                                                            #    13.0   703.8  1233.5  1379.9  1963.5  4147.0

    save(GTEX_c, file=paste0(output_dir,"GTEX_c",".Rdata"))


}else if(opt$database=="CCLE"){
    ensg_upe1<-paste(ensg_upe1_ccle)
    ensg_pe1_sin_upe1<-paste(ensg_pe1_sin_upe1_ccle)

    CCLE_sel_upe1<-intersect(ensg_upe1, rownames(ccle_filtrado_5_25_norm))
    CCLE_sel_PE1<-intersect(ensg_pe1_sin_upe1, rownames(ccle_filtrado_5_25_norm))

    CCLE_cor <- apply(ccle_filtrado_5_25_norm[CCLE_sel_upe1,], 1, FUN=function(x) apply(ccle_filtrado_5_25_norm[CCLE_sel_PE1,], 1, FUN=function(y) cor(x,y)))
    write.table(rownames(CCLE_cor), file=paste0(output_dir,"CCLE_PE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")
    write.table(colnames(CCLE_cor), file=paste0(output_dir,"CCLE_UPE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")

    save(CCLE_cor, file=paste0(output_dir,"CCLE_cor",".Rdata"))
    CCLE_corP <- apply(ccle_filtrado_5_25_norm[CCLE_sel_upe1,], 1, FUN=function(x) apply(ccle_filtrado_5_25_norm[CCLE_sel_PE1,], 1, FUN=function(y) cor.test(x,y)$p.value))
    save(CCLE_corP, file=paste0(output_dir,"CCLE_corP",".Rdata"))
    CCLE_corFDR <- apply(CCLE_corP, 2, FUN=function(x) p.adjust(x,method="fdr"))
    save(CCLE_corFDR, file=paste0(output_dir,"CCLE_corFDR",".Rdata"))
    CCLE_c <- (apply(CCLE_corFDR, 2, FUN=function(x) (x<0.05)*1)) * (apply(CCLE_cor, 2, FUN=function(x) (abs(x)>0.3)*1)) # meanGS=555.1769, medianGS=233.5, 0_GS=25 #eli
                                                                                                                         #   summary(colSums(data))
                                                                                                                         #      Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
                                                                                                                         #       0.0    36.0   197.0   512.0   725.2  4060.0   0_GS=55      

    save(CCLE_c, file=paste0(output_dir,"CCLE_c",".Rdata"))

                                                                                                                        ###PE 10 % 
                                                                                                                        #summary(colSums(CCLE_c))
                                                                                                                        #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
                                                                                                                        #1.0   167.0   501.5   752.9  1062.0  3781.0 
                                                                                                                        #sum(colSums(CCLE_c)==0)
                                                                                                                        #0
}
name<-basename(opt$file)
name_s<-unlist(lapply(strsplit(paste(name),"_"),"[",2))


mkdir(paste0(output_dir,"chunks_",name_s,"/"))
output_chunks<-paste0(output_dir,"/chunks_",name_s,"/")

if(name_s =="GOBP"){   #separar por columnas
    geneXGO<-as.data.frame(geneXGO)
    n <- 300
    nc <- ncol(geneXGO)
    geneXGO_splitted<-list()
    k<-1
    max_len<-length(names(geneXGO))
    for(i in seq(from=1, to=max_len, by=n)){


        # last loop
        if((i)+n > (max_len)){
        chunk<-geneXGO[,c((i):max_len)]
        geneXGO_splitted[[k]]<-chunk
        cat(i,"_",max_len,"\n")
    }
        else{
        chunk<-geneXGO[,c((i):(i+n-1))]
        geneXGO_splitted[[k]]<-chunk
        k<-k+1
        cat((i),"_",(i+n-1),"\n")
    }

    }
    for(i in 1:length(geneXGO_splitted)){

        cat(i,"\n")
        chunk<-geneXGO_splitted[[i]]
        save(chunk, file =paste0(output_chunks,"geneXGO_chunk_",i,".Rdata"))

    }

}else if(name_s=="GOCC"){

    geneXGO_cc<-as.data.frame(geneXGO_CC)
    n <- 300
    nc <- ncol(geneXGO_cc)
    geneXGO_cc_splitted<-list()
    k<-1
    max_len<-length(names(geneXGO_cc))
    for(i in seq(from=1, to=max_len, by=n)){
        cat(i,"\n")
        # last loop
        if(i+n > max_len){
        chunk<-geneXGO_cc[,c((i):max_len)]
        geneXGO_cc_splitted[[k]]<-chunk
    }
        else{
        chunk<-geneXGO_cc[,c((i):(i+n-1))]
        geneXGO_cc_splitted[[k]]<-chunk
        k<-k+1

        }
    }
    for(i in 1:length(geneXGO_cc_splitted)){

        cat(i,"\n")
        chunk<-geneXGO_cc_splitted[[i]]
        save(chunk, file =paste0(output_chunks,"geneXGO_CC_chunk_",i,".Rdata"))

    }

}else if(name_s=="GOMF"){

    geneXGO_MF<-as.data.frame(geneXGO_MF)
    n <- 300
    nc <- ncol(geneXGO_MF)
    geneXGO_MF_splitted<-list()
    k<-1
    max_len<-length(names(geneXGO_MF))
    for(i in seq(from=1, to=max_len, by=n)){
        cat(i,"\n")
        # last loop
        if(i+n > max_len){
        chunk<-geneXGO_MF[,c((i):max_len)]
        geneXGO_MF_splitted[[k]]<-chunk
    }
        else{
        chunk<-geneXGO_MF[,c((i):(i+n-1))]
        geneXGO_MF_splitted[[k]]<-chunk
        k<-k+1

        }
    }
    for(i in 1:length(geneXGO_MF_splitted)){

        cat(i,"\n")
        chunk<-geneXGO_MF_splitted[[i]]
        save(chunk, file =paste0(output_chunks,"geneXGO_MF_chunk_",i,".Rdata"))

    }

}else if(name_s=="MSigDB"){

    msigdb_mat_ensg2<-as.data.frame(msigdb_mat_ensg2)
    n <- 300
    nc <- ncol(msigdb_mat_ensg2)
    msigdb_mat_ensg2_splitted<-list()
    k<-1
    max_len<-length(names(msigdb_mat_ensg2))
    for(i in seq(from=1, to=max_len, by=n)){
        cat(i,"\n")
        # last loop
        if(i+n > max_len){
        chunk<-msigdb_mat_ensg2[,c((i):max_len)]
        msigdb_mat_ensg2_splitted[[k]]<-chunk
    }
        else{
        chunk<-msigdb_mat_ensg2[,c((i):(i+n-1))]
        msigdb_mat_ensg2_splitted[[k]]<-chunk
        k<-k+1

        }
    }
    for(i in 1:length(msigdb_mat_ensg2_splitted)){

        cat(i,"\n")
        chunk<-msigdb_mat_ensg2_splitted[[i]]
        save(chunk, file =paste0(output_chunks,"msigdb_mat_ensg2_chunk_",i,".Rdata"))

    }

}else if(name_s=="Disease"){

    #malaCards_mat_ensg2<-as.data.frame(malaCards_mat_ensg2)
    #colnames(malaCards_mat_ensg2)<-gsub(" ","@", colnames(malaCards_mat_ensg2))
    #malaCards_mat_ensg2<-as.data.frame(malaCards_mat_ensg2)
    n <- 300
    nc <- ncol(malaCards_mat_ensg2)
    malaCards_mat_ensg2_splitted<-list()
    k<-1
    max_len<-length(colnames(malaCards_mat_ensg2))
    for(i in seq(from=1, to=max_len, by=n)){
        cat(i,"\n")
        # last loop
        if(i+n > max_len){
        chunk<-malaCards_mat_ensg2[,c((i):max_len)]
        malaCards_mat_ensg2_splitted[[k]]<-chunk
        k<-k+1
    }
        else{
        chunk<-malaCards_mat_ensg2[,c((i):(i+n-1))]
        malaCards_mat_ensg2_splitted[[k]]<-chunk
        k<-k+1

        }
    }
    for(i in 1:length(malaCards_mat_ensg2_splitted)){

        cat(i,"\n")
        chunk<-malaCards_mat_ensg2_splitted[[i]]
        save(chunk, file =paste0(output_chunks,"malaCards_mat_ensg2_chunk_",i,".Rdata"))

    }


}



cat("done step1,\n")
