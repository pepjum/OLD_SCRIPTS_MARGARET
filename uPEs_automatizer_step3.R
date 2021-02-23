### step 3
args=(commandArgs(TRUE))

folder_output<-args[1] #OUTPUT_DATOS_DATABASE
datos<-args[2] ## CCLE, TCGA ....
database<-args[3]  #GOBP, malaCards, GO_cc


files<-list.files(folder_output, pattern=".Rdata")
files<-paste0(folder_output,files)

folder_less<-dirname(folder_output)
cat("folder_less ", folder_less,"\n")

if(datos=="TCGA"){

    TCGA_uPE1_enrichedGO_pval<-data.frame()
    TCGA_uPE1_enrichedGO_geneXenrichedGO<-data.frame()
    for(i in 1:length(files)){
        cat(i,"\n")
        loaded<-get(load(files[[i]]))

        TCGA_uPE1_enrichedGO_pval<-rbind(TCGA_uPE1_enrichedGO_pval,loaded$pval)
        TCGA_uPE1_enrichedGO_geneXenrichedGO<-rbind(TCGA_uPE1_enrichedGO_geneXenrichedGO, loaded$geneXenrichedGO)
        cat(dim(TCGA_uPE1_enrichedGO_pval),"\n")
    }
    TCGA_uPE1_enrichedGO_pval_sel2<-TCGA_uPE1_enrichedGO_pval

    TCGA_uPE1_enrichedGO_pval_sel <- (apply(TCGA_uPE1_enrichedGO_pval,2,FUN=function(x) (x<0.01)*1))
    TCGA_uPE1_enrichedGO_pval_sel[is.na(TCGA_uPE1_enrichedGO_pval_sel)] <- 0
   
    
    if(database=="GOBP"|| database=="GOCC" || database=="GOMF"){
        validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/GOannotation_Ago20.txt", header=T, as.is=T, quote="", sep="\t", fill=T)
        TCGA_uPE1_enrichedGO_pval_sel_f<-TCGA_uPE1_enrichedGO_pval_sel[rownames(TCGA_uPE1_enrichedGO_pval_sel) %in% validated$GOID,]
        TCGA_uPE1_enrichedGO_pval_sel_f2<-TCGA_uPE1_enrichedGO_pval_sel2[rownames(TCGA_uPE1_enrichedGO_pval_sel2) %in% validated$GOID,]
        
        TCGA_uPE1_enrichedGO_geneXenrichedGO_f<-TCGA_uPE1_enrichedGO_geneXenrichedGO[TCGA_uPE1_enrichedGO_geneXenrichedGO$GO %in% validated$GOID,]
        
    }else if(database=="Disease"){
        # validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/malaCards_annotation.txt", header=T, as.is=T, quote="", sep="\t")
        TCGA_uPE1_enrichedGO_pval_sel_f<-TCGA_uPE1_enrichedGO_pval_sel
        TCGA_uPE1_enrichedGO_pval_sel_f2<-TCGA_uPE1_enrichedGO_pval_sel2
        
        TCGA_uPE1_enrichedGO_geneXenrichedGO_f<-TCGA_uPE1_enrichedGO_geneXenrichedGO

    }else if(database=="MSigDB"){
        # validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/msigdb_annotation.txt", header=T, as.is=T, quote="", sep="\t")
        TCGA_uPE1_enrichedGO_pval_sel_f<-TCGA_uPE1_enrichedGO_pval_sel
        TCGA_uPE1_enrichedGO_pval_sel_f2<-TCGA_uPE1_enrichedGO_pval_sel2
        
        TCGA_uPE1_enrichedGO_geneXenrichedGO_f<-TCGA_uPE1_enrichedGO_geneXenrichedGO

    }
    TCGA_uPE1_enrichedGO_geneXenrichedGO<-TCGA_uPE1_enrichedGO_geneXenrichedGO_f

    write.table(cbind(GO_id=rownames(TCGA_uPE1_enrichedGO_pval_sel_f),TCGA_uPE1_enrichedGO_pval_sel_f), file=paste0(folder_output,"TCGA_enriched_",database,".txt"), row.names=F, sep="\t", quote=F)
    write.table(cbind(GO_id=rownames(TCGA_uPE1_enrichedGO_pval_sel_f2),TCGA_uPE1_enrichedGO_pval_sel_f2), file=paste0(folder_output,"TCGA_enriched_PVAL_",database,".txt"), row.names=F, sep="\t", quote=T)

    write.table(rownames(TCGA_uPE1_enrichedGO_pval_sel_f), file=paste0(folder_output,"TCGA_",database,"_names.txt"), row.names=F, col.names=F, sep="\t", quote=F)

    TCGA_c<-get(load(paste0(folder_less,"/","TCGA_c.Rdata")))

    TCGA_c_order <- TCGA_c[intersect(rownames(TCGA_c),TCGA_uPE1_enrichedGO_geneXenrichedGO[,1]),]
    TCGA_c_order <- TCGA_c_order[order(rownames(TCGA_c_order)),]
    TCGA_c_order <- rbind(TCGA_c_order,TCGA_c[setdiff(rownames(TCGA_c),TCGA_uPE1_enrichedGO_geneXenrichedGO[,1]),])

    # TCGA_GOCC2 <-TCGA_GOCC[,-1]
    # rownames(TCGA_GOCC2) <- paste(TCGA_GOCC[,1])
    # TCGA_GOBP2 <- TCGA_GOBP[,-1]
    # rownames(TCGA_GOBP2) <- TCGA_GOBP[,1]
    #
    # TCGA_GOCC2_o <- TCGA_GOCC2[rownames(TCGA_GOBP2),colnames(TCGA_GOBP2)]


# > TCGA_GeneXenriched2_ord<-TCGA_GeneXenriched2[rownames(TCGA_GeneXenriched_GOBP2)]
# Error in `[.data.frame`(TCGA_GeneXenriched2, rownames(TCGA_GeneXenriched_GOBP2)) :
#   undefined columns selected
# > TCGA_GeneXenriched2_ord<-TCGA_GeneXenriched2[rownames(TCGA_GeneXenriched_GOBP2),]
# > dim(TCGA_GeneXenriched2_ord)
# [1] 15288  4985
# > dim(TCGA_GeneXenriched_GOBP)
# [1] 15288 11855
# > sum(rownames(TCGA_GeneXenriched2_ord) == rownames(TCGA_GeneXenriched_GOBP2))
# [1]


    write.table(cbind(PE1=rownames(TCGA_c_order),TCGA_c_order), file=paste0(folder_output,"TCGA_",database,"_corMatrix.txt"), row.names=F, sep="\t", quote=F)   # cuando genere este, usarlo como referencia para el resto
    # write.table(TCGA_uPE1_enrichedGO_pval_sel_f$GO_id, file=paste0(folder_output,"TCGA_",database,"_names.txt"), row.names=F, col.names=F, sep="\n", quote=F)
    write.table(rownames(TCGA_c_order), file=paste0(folder_output,"TCGA_",database,"_PE_names.txt"), row.names=F, col.names=F, quote=F, sep="\n")
    write.table(colnames(TCGA_c_order), file=paste0(folder_output,"TCGA_",database,"_UPE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")


        ##### por aqui

    TCGA_geneXgo_order <- TCGA_uPE1_enrichedGO_geneXenrichedGO[TCGA_uPE1_enrichedGO_geneXenrichedGO[,1] %in% intersect(rownames(TCGA_c),TCGA_uPE1_enrichedGO_geneXenrichedGO[,1]),]
    TCGA_geneXgo_order <- TCGA_geneXgo_order[order(TCGA_geneXgo_order[,1]),]
    library(reshape2)
    TCGA_geneXgo_order_mat <- dcast(TCGA_geneXgo_order, Gene~GO)
    TCGA_geneXgo_order_mat <- TCGA_geneXgo_order_mat[order(paste(TCGA_geneXgo_order_mat[,1])),]
    TCGA_geneWoGO <- t(sapply(setdiff(rownames(TCGA_c),TCGA_uPE1_enrichedGO_geneXenrichedGO[,1]), FUN=function(x) c(x, rep(NA,ncol(TCGA_geneXgo_order_mat)-1))))
    colnames(TCGA_geneWoGO) <- colnames(TCGA_geneXgo_order_mat)
    TCGA_geneXgo_order_mat <- rbind(TCGA_geneXgo_order_mat, TCGA_geneWoGO)
    TCGA_geneXgo_order_mat01 <- apply(TCGA_geneXgo_order_mat[,-1], 2, FUN=function(x) ifelse(is.na(x),0,1))
    rownames(TCGA_geneXgo_order_mat01) <- paste(TCGA_geneXgo_order_mat[,1])
    TCGA_GOWoGene <-  sapply(setdiff(rownames(TCGA_uPE1_enrichedGO_pval_sel_f),colnames(TCGA_geneXgo_order_mat01)), FUN=function(x) rep(0,nrow(TCGA_geneXgo_order_mat01)))
    rownames(TCGA_GOWoGene) <- rownames(TCGA_geneXgo_order_mat01)
    TCGA_geneXgo_order_mat01 <- cbind(TCGA_geneXgo_order_mat01, TCGA_GOWoGene)

    TCGA_geneXgo_order_mat01 <- TCGA_geneXgo_order_mat01[,rownames(TCGA_uPE1_enrichedGO_pval_sel_f)]
    write.table(cbind(Gene=rownames(TCGA_geneXgo_order_mat01),TCGA_geneXgo_order_mat01), file=paste0(folder_output,"TCGA_GeneXenriched_",database,".txt"), row.names=F, sep="\t", quote=F)



}else if (datos=="GTEX"){

    GTEX_uPE1_enrichedGO_pval<-data.frame()
    GTEX_uPE1_enrichedGO_geneXenrichedGO<-data.frame()
    for(i in 1:length(files)){
        cat(i,"\n")
        loaded<-get(load(files[[i]]))

        GTEX_uPE1_enrichedGO_pval<-rbind(GTEX_uPE1_enrichedGO_pval,loaded$pval)
        GTEX_uPE1_enrichedGO_geneXenrichedGO<-rbind(GTEX_uPE1_enrichedGO_geneXenrichedGO, loaded$geneXenrichedGO)
        cat(dim(GTEX_uPE1_enrichedGO_pval),"\n")
    }

    GTEX_uPE1_enrichedGO_pval_sel2<-GTEX_uPE1_enrichedGO_pval
    GTEX_uPE1_enrichedGO_pval_sel <- (apply(GTEX_uPE1_enrichedGO_pval,2,FUN=function(x) (x<0.01)*1))
    GTEX_uPE1_enrichedGO_pval_sel[is.na(GTEX_uPE1_enrichedGO_pval_sel)] <- 0

    if(database=="GOBP"|| database=="GOCC" || database=="GOMF"){
        validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/GOannotation_Ago20.txt", header=T, as.is=T, quote="", sep="\t", fill=T)
        GTEX_uPE1_enrichedGO_pval_sel_f<-GTEX_uPE1_enrichedGO_pval_sel[rownames(GTEX_uPE1_enrichedGO_pval_sel) %in% validated$GOID,]
        GTEX_uPE1_enrichedGO_pval_sel_f2<-GTEX_uPE1_enrichedGO_pval_sel2[rownames(GTEX_uPE1_enrichedGO_pval_sel2) %in% validated$GOID,]
        

        GTEX_uPE1_enrichedGO_geneXenrichedGO_f<-GTEX_uPE1_enrichedGO_geneXenrichedGO[GTEX_uPE1_enrichedGO_geneXenrichedGO$GO %in% validated$GOID,]

    }else if(database=="Disease"){
        # validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/malaCards_annotation.txt", header=T, as.is=T, quote="", sep="\t")
        GTEX_uPE1_enrichedGO_pval_sel_f<-GTEX_uPE1_enrichedGO_pval_sel
        GTEX_uPE1_enrichedGO_pval_sel_f2<-GTEX_uPE1_enrichedGO_pval_sel2
        
        GTEX_uPE1_enrichedGO_geneXenrichedGO_f<-GTEX_uPE1_enrichedGO_geneXenrichedGO

    }else if(database=="MSigDB"){
        # validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/msigdb_annotation.txt", header=T, as.is=T, quote="", sep="\t")
        GTEX_uPE1_enrichedGO_pval_sel_f<-GTEX_uPE1_enrichedGO_pval_sel
        GTEX_uPE1_enrichedGO_pval_sel_f2<-GTEX_uPE1_enrichedGO_pval_sel2
        
        GTEX_uPE1_enrichedGO_geneXenrichedGO_f<-GTEX_uPE1_enrichedGO_geneXenrichedGO

    }
    GTEX_uPE1_enrichedGO_geneXenrichedGO<-GTEX_uPE1_enrichedGO_geneXenrichedGO_f

    write.table(cbind(GO_id=rownames(GTEX_uPE1_enrichedGO_pval_sel_f),GTEX_uPE1_enrichedGO_pval_sel_f), file=paste0(folder_output,"GTEX_enriched_",database,".txt"), row.names=F, sep="\t", quote=F)
    write.table(cbind(GO_id=rownames(GTEX_uPE1_enrichedGO_pval_sel_f2),GTEX_uPE1_enrichedGO_pval_sel_f2), file=paste0(folder_output,"GTEX_enriched_",database,"_PVAL.txt"), row.names=F, sep="\t", quote=T)


    write.table(rownames(GTEX_uPE1_enrichedGO_pval_sel_f), file=paste0(folder_output,"GTEX_",database,"_names.txt"), row.names=F, col.names=F, sep="\n", quote=F)

    GTEX_c<-get(load(paste0(folder_less,"/","GTEX_c.Rdata")))

    GTEX_c_order <- GTEX_c[intersect(rownames(GTEX_c),GTEX_uPE1_enrichedGO_geneXenrichedGO[,1]),]
    GTEX_c_order <- GTEX_c_order[order(rownames(GTEX_c_order)),]
    GTEX_c_order <- rbind(GTEX_c_order,GTEX_c[setdiff(rownames(GTEX_c),GTEX_uPE1_enrichedGO_geneXenrichedGO[,1]),])

    write.table(cbind(PE1=rownames(GTEX_c_order),GTEX_c_order), file=paste0(folder_output,"GTEX_",database,"_corMatrix.txt"), row.names=F, sep="\t", quote=F)
    # //write.table(GTEX_uPE1_enrichedGO_pval_sel_f$GO_id, file=paste0(folder_output,"GTEX_",database,"_names.txt"), row.names=F, col.names=F, sep="\n", quote=F)

     # write.table(rownames(GTEX_c_order), file=paste0(folder_output,"GTEX_",database,"_names.txt"), row.names=F, col.names=F, sep="\n", quote=F)
    write.table(rownames(GTEX_c_order), file=paste0(folder_output,"GTEX_",database,"_PE_names.txt"), row.names=F, col.names=F, quote=F, sep="\n")
    write.table(colnames(GTEX_c_order), file=paste0(folder_output,"GTEX_",database,"_UPE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")

        ##### por aqui

    GTEX_geneXgo_order <- GTEX_uPE1_enrichedGO_geneXenrichedGO[GTEX_uPE1_enrichedGO_geneXenrichedGO[,1] %in% intersect(rownames(GTEX_c),GTEX_uPE1_enrichedGO_geneXenrichedGO[,1]),]
    GTEX_geneXgo_order <- GTEX_geneXgo_order[order(GTEX_geneXgo_order[,1]),]
    library(reshape2)
    GTEX_geneXgo_order_mat <- dcast(GTEX_geneXgo_order, Gene~GO)
    GTEX_geneXgo_order_mat <- GTEX_geneXgo_order_mat[order(paste(GTEX_geneXgo_order_mat[,1])),]
    GTEX_geneWoGO <- t(sapply(setdiff(rownames(GTEX_c),GTEX_uPE1_enrichedGO_geneXenrichedGO[,1]), FUN=function(x) c(x, rep(NA,ncol(GTEX_geneXgo_order_mat)-1))))
    colnames(GTEX_geneWoGO) <- colnames(GTEX_geneXgo_order_mat)
    GTEX_geneXgo_order_mat <- rbind(GTEX_geneXgo_order_mat, GTEX_geneWoGO)
    GTEX_geneXgo_order_mat01 <- apply(GTEX_geneXgo_order_mat[,-1], 2, FUN=function(x) ifelse(is.na(x),0,1))
    rownames(GTEX_geneXgo_order_mat01) <- paste(GTEX_geneXgo_order_mat[,1])
    GTEX_GOWoGene <-  sapply(setdiff(rownames(GTEX_uPE1_enrichedGO_pval_sel_f),colnames(GTEX_geneXgo_order_mat01)), FUN=function(x) rep(0,nrow(GTEX_geneXgo_order_mat01)))
    rownames(GTEX_GOWoGene) <- rownames(GTEX_geneXgo_order_mat01)
    GTEX_geneXgo_order_mat01 <- cbind(GTEX_geneXgo_order_mat01, GTEX_GOWoGene)

    GTEX_geneXgo_order_mat01 <- GTEX_geneXgo_order_mat01[,rownames(GTEX_uPE1_enrichedGO_pval_sel_f)]
    write.table(cbind(Gene=rownames(GTEX_geneXgo_order_mat01),GTEX_geneXgo_order_mat01), file=paste0(folder_output,"GTEX_GeneXenriched_",database,".txt"), row.names=F, sep="\t", quote=F)


}else if(datos=="CCLE"){

    CCLE_uPE1_enrichedGO_pval<-data.frame()
    CCLE_uPE1_enrichedGO_geneXenrichedGO<-data.frame()
    for(i in 1:length(files)){
        cat(i,"\n")
        loaded<-get(load(files[[i]]))

        CCLE_uPE1_enrichedGO_pval<-rbind(CCLE_uPE1_enrichedGO_pval,loaded$pval)
        
        CCLE_uPE1_enrichedGO_geneXenrichedGO<-rbind(CCLE_uPE1_enrichedGO_geneXenrichedGO, loaded$geneXenrichedGO)
        cat(dim(CCLE_uPE1_enrichedGO_pval),"\n")
    }
    CCLE_uPE1_enrichedGO_pval_sel2 <- CCLE_uPE1_enrichedGO_pval
    CCLE_uPE1_enrichedGO_pval_sel <- (apply(CCLE_uPE1_enrichedGO_pval,2,FUN=function(x) (x<0.01)*1))
    CCLE_uPE1_enrichedGO_pval_sel[is.na(CCLE_uPE1_enrichedGO_pval_sel)] <- 0
    
    if(database=="GOBP"|| database=="GOCC" || database=="GOMF"){
        validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/GOannotation_Ago20.txt", header=T, as.is=T, quote="", sep="\t", fill=T)
        CCLE_uPE1_enrichedGO_pval_sel_f<-CCLE_uPE1_enrichedGO_pval_sel[rownames(CCLE_uPE1_enrichedGO_pval_sel) %in% validated$GOID,]
        CCLE_uPE1_enrichedGO_pval_sel_f2<-CCLE_uPE1_enrichedGO_pval_sel2[rownames(CCLE_uPE1_enrichedGO_pval_sel2) %in% validated$GOID,]
 
        CCLE_uPE1_enrichedGO_geneXenrichedGO_f<-CCLE_uPE1_enrichedGO_geneXenrichedGO[CCLE_uPE1_enrichedGO_geneXenrichedGO$GO %in% validated$GOID,]

    }else if(database=="Disease"){
#        validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/malaCards_annotation.txt", header=T, as.is=T, quote="", sep="\t")
        CCLE_uPE1_enrichedGO_pval_sel_f<-CCLE_uPE1_enrichedGO_pval_sel
        CCLE_uPE1_enrichedGO_pval_sel_f2<-CCLE_uPE1_enrichedGO_pval_sel2
        
        CCLE_uPE1_enrichedGO_geneXenrichedGO_f<-CCLE_uPE1_enrichedGO_geneXenrichedGO

    }else if(database=="MSigDB"){
#        validated<-read.table("~/data/pepe/59_uPEs_feb20/GuillermoDB/msigdb_annotation.txt", header=T, as.is=T, quote="", sep="\t")
        CCLE_uPE1_enrichedGO_pval_sel_f<-CCLE_uPE1_enrichedGO_pval_sel
        CCLE_uPE1_enrichedGO_pval_sel_f2<-CCLE_uPE1_enrichedGO_pval_sel2
        
        CCLE_uPE1_enrichedGO_geneXenrichedGO_f<-CCLE_uPE1_enrichedGO_geneXenrichedGO

    }
    CCLE_uPE1_enrichedGO_geneXenrichedGO<-CCLE_uPE1_enrichedGO_geneXenrichedGO_f
    cat("outputs","\n")
    write.table(cbind(GO_id=rownames(CCLE_uPE1_enrichedGO_pval_sel_f),CCLE_uPE1_enrichedGO_pval_sel_f), file=paste0(folder_output,"CCLE_enriched_",database,".txt"), row.names=F, sep="\t", quote=F)
    write.table(cbind(GO_id=rownames(CCLE_uPE1_enrichedGO_pval_sel_f2),CCLE_uPE1_enrichedGO_pval_sel_f2), file=paste0(folder_output,"CCLE_enriched_",database,"_PVAL.txt"), row.names=F, sep="\t", quote=T, col.names=T)
    
    write.table(rownames(CCLE_uPE1_enrichedGO_pval_sel_f), file=paste0(folder_output,"CCLE_",database,"_names.txt"), row.names=F, col.names=F, sep="\n", quote=F)

    CCLE_c<-get(load(paste0(folder_less,"/","CCLE_c.Rdata")))

    CCLE_c_order <- CCLE_c[intersect(rownames(CCLE_c),CCLE_uPE1_enrichedGO_geneXenrichedGO[,1]),]
    CCLE_c_order <- CCLE_c_order[order(rownames(CCLE_c_order)),]
    CCLE_c_order <- rbind(CCLE_c_order,CCLE_c[setdiff(rownames(CCLE_c),CCLE_uPE1_enrichedGO_geneXenrichedGO[,1]),])

    write.table(cbind(PE1=rownames(CCLE_c_order),CCLE_c_order), file=paste0(folder_output,"CCLE_",database,"_corMatrix.txt"), row.names=F, sep="\t", quote=F)
#    write.table(rownames(CCLE_c_order), file=paste0(folder_output,"CCLE_",database,"_names.txt"), row.names=F, col.names=F, sep="\n", quote=F)
    # write.table(CCLE_uPE1_enrichedGO_pval_sel_f$GO_id, file=paste0(folder_output,"CCLE_",database,"_names.txt"), row.names=F, col.names=F, sep="\n", quote=F)


    write.table(rownames(CCLE_c_order), file=paste0(folder_output,"CCLE_",database,"_PE_names.txt"), row.names=F, col.names=F, quote=F, sep="\n")
    write.table(colnames(CCLE_c_order), file=paste0(folder_output,"CCLE_",database,"_UPE_names.txt"),row.names=F, col.names=F, quote=F, sep="\n")

        ##### por aqui

    CCLE_geneXgo_order <- CCLE_uPE1_enrichedGO_geneXenrichedGO[CCLE_uPE1_enrichedGO_geneXenrichedGO[,1] %in% intersect(rownames(CCLE_c),CCLE_uPE1_enrichedGO_geneXenrichedGO[,1]),]
    CCLE_geneXgo_order <- CCLE_geneXgo_order[order(CCLE_geneXgo_order[,1]),]
    library(reshape2)
    CCLE_geneXgo_order_mat <- dcast(CCLE_geneXgo_order, Gene~GO)
    CCLE_geneXgo_order_mat <- CCLE_geneXgo_order_mat[order(paste(CCLE_geneXgo_order_mat[,1])),]
    CCLE_geneWoGO <- t(sapply(setdiff(rownames(CCLE_c),CCLE_uPE1_enrichedGO_geneXenrichedGO[,1]), FUN=function(x) c(x, rep(NA,ncol(CCLE_geneXgo_order_mat)-1))))
    colnames(CCLE_geneWoGO) <- colnames(CCLE_geneXgo_order_mat)
    CCLE_geneXgo_order_mat <- rbind(CCLE_geneXgo_order_mat, CCLE_geneWoGO)
    CCLE_geneXgo_order_mat01 <- apply(CCLE_geneXgo_order_mat[,-1], 2, FUN=function(x) ifelse(is.na(x),0,1))
    rownames(CCLE_geneXgo_order_mat01) <- paste(CCLE_geneXgo_order_mat[,1])
    CCLE_GOWoGene <-  sapply(setdiff(rownames(CCLE_uPE1_enrichedGO_pval_sel_f),colnames(CCLE_geneXgo_order_mat01)), FUN=function(x) rep(0,nrow(CCLE_geneXgo_order_mat01)))
    rownames(CCLE_GOWoGene) <- rownames(CCLE_geneXgo_order_mat01)
    CCLE_geneXgo_order_mat01 <- cbind(CCLE_geneXgo_order_mat01, CCLE_GOWoGene)

    CCLE_geneXgo_order_mat01 <- CCLE_geneXgo_order_mat01[,rownames(CCLE_uPE1_enrichedGO_pval_sel_f)]
    write.table(cbind(Gene=rownames(CCLE_geneXgo_order_mat01),CCLE_geneXgo_order_mat01), file=paste0(folder_output,"CCLE_GeneXenriched_",database,".txt"), row.names=F, sep="\t", quote=F)



}
cat("DONE!","\n")
#### El orden de las filas de los ficheros geneXenriched de GTEX debe ser igual entre si. Igual para CCLE y para TCGA_
#### El orden de las columnas de los ficheros enriched de GTEX debe ser igual entre si. Lo mismo para las otras bases de datos
#### Se van a generar 5 ficheros cormatrix con distintos nombres ... GOBP, GOCC, malacards, etc... pero su contenido es el mismo. Solo se necesita 1
### los ficheros PE names y UPE names se van a generar igual que los cormatrix, y deben coincidir en orden con la columna o fila correspondiente de los geneenriched o geneXenrichedGO y con los _c.Rdata
