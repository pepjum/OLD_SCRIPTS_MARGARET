

#############################################################################################
################################################################################################
NUEVOS ORDENAMIENTOS

##################################
#### TCGA

### check GOBP
library(tidyverse)

TCGA_enriched_GOBP<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_enriched_GOBP.txt", col_names=T)

TCGA_geneXenriched_GOBP<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_GeneXenriched_GOBP.txt", col_names=T)
#sum(colnames(TCGA_geneXenriched_GOBP2) == TCGA_enriched_GOBP$GO_id)
#[1] 11854
#> dim(TCGA_geneXenriched_GOBP2)
#[1] 15296 11854
#> dim(TCGA_enriched_GOBP2)
#[1] 11854  1297

cormatrix<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_GOBP_corMatrix.txt")
#> sum(cormatrix$PE1 == TCGA_geneXenriched_GOBP$Gene)

##### GOCC ###############################################3
### ahora hay que ordenar GOCC en funcion de los names que ha generado  GOBP

TCGA_enriched_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_enriched_GOCC.txt", col_names=T)

TCGA_geneXenriched_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_GeneXenriched_GOCC.txt", col_names=T)

#TCGA_enriched_GOCC2<-TCGA_enriched_GOCC[,-1]
# dim(TCGA_enriched_GOCC2)
#[1] 1339 1297
# head(TCGA_enriched_GOCC2[1:2,1:2])
# A tibble: 2 x 2
#  ENSG00000000460 ENSG00000003249
#            <dbl>           <dbl>
#1               0               0
#2               0               0
# sum(colnames(TCGA_enriched_GOCC2) == colnames(TCGA_enriched_GOBP2))
#[1] 1297
# sum(TCGA_geneXenriched_GOBP$Gene == TCGA_geneXenriched_GOCC$Gene)
### 172 (hay que reordenar)
TCGA_geneXenriched_GOCC2<-TCGA_geneXenriched_GOCC[order(match(TCGA_geneXenriched_GOCC$Gene, TCGA_geneXenriched_GOBP$Gene)),]

write_tsv(TCGA_geneXenriched_GOCC2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_GeneXenriched_GOCC.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de GOCC en funcion de TCGA_geneXenriched_GOCC$Gene

cormatrix_TCGA_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_GOCC_corMatrix.txt", col_names=T)
#sum(cormatrix_TCGA_GOCC$PE1 == TCGA_geneXenriched_GOCC2$Gene)
#172
cormatrix_TCGA_GOCC2<-cormatrix_TCGA_GOCC[order(match(cormatrix_TCGA_GOCC$PE1, TCGA_geneXenriched_GOCC2$Gene)),]
#sum(cormatrix_TCGA_GOCC2$PE1 == TCGA_geneXenriched_GOCC2$Gene)
write.table(cormatrix_TCGA_GOCC2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_GOCC_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(cormatrix_TCGA_GOCC2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")


#comprobar los uPE (colnames)
cormatrix_TCGA_GOCC3<-cormatrix_TCGA_GOCC2[,-1]
#sum(colnames(TCGA_enriched_GOBP2 == colnames(cormatrix_TCGA_GOCC3)))
#1297
write.table(colnames(cormatrix_TCGA_GOCC3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_GOCC_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(colnames(cormatrix_TCGA_GOCC3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

##### ORDENAR TCGA_C y TCGA_corFDR para que tambien tengan el mismo orden

TCGA_c<-get(load("~/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_c.Rdata"))
TCGA_c2<-TCGA_c[order(match(rownames(TCGA_c), cormatrix_TCGA_GOCC2$PE1)),]
sum(colnames(TCGA_enriched_GOBP2)== colnames(TCGA_c2))
save(TCGA_c2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_c.Rdata")

TCGA_cor<-get(load("~/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_cor.Rdata"))
TCGA_cor2<-TCGA_cor[order(match(rownames(TCGA_cor), cormatrix_TCGA_GOCC2$PE1)),]
save(TCGA_cor2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_cor.Rdata")

TCGA_corP<-get(load("~/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_corP.Rdata"))
TCGA_corP2<-TCGA_corP[order(match(rownames(TCGA_corP), cormatrix_TCGA_GOCC2$PE1)),]
save(TCGA_corP2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_corP.Rdata")

TCGA_corFDR<-get(load("~/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_corFDR.Rdata"))
TCGA_corFDR2<-TCGA_corFDR[order(match(rownames(TCGA_corFDR), cormatrix_TCGA_GOCC2$PE1)),]
save(TCGA_corFDR2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/TCGA_corFDR.Rdata")

####### ORDENAR GOMF

#TCGA_enriched_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_enriched_GOMF.txt", col_names=T)

TCGA_geneXenriched_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_GeneXenriched_GOMF.txt", col_names=T)

#TCGA_enriched_GOMF2<-TCGA_enriched_GOMF[,-1]

TCGA_geneXenriched_GOMF2<-TCGA_geneXenriched_GOMF[order(match(TCGA_geneXenriched_GOMF$Gene, TCGA_geneXenriched_GOBP$Gene)),]

write_tsv(TCGA_geneXenriched_GOMF2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_GeneXenriched_GOMF.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de GOMF en funcion de TCGA_geneXenriched_GOMF$Gene

cormatrix_TCGA_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_GOMF_corMatrix.txt")

cormatrix_TCGA_GOMF2<-cormatrix_TCGA_GOMF[order(match(cormatrix_TCGA_GOMF$PE1, TCGA_geneXenriched_GOMF2$Gene)),]

write.table(cormatrix_TCGA_GOMF2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_GOMF_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_TCGA_GOMF3<-cormatrix_TCGA_GOMF2[,-1]
#sum(colnames(TCGA_enriched_GOBP2) == colnames(cormatrix_TCGA_GOMF3))
#1297
write.table(colnames(cormatrix_TCGA_GOMF3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_GOMF_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

####### ORDENAR Disease

TCGA_enriched_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_enriched_Disease.txt", col_names=T)

TCGA_geneXenriched_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_GeneXenriched_Disease.txt", col_names=T)

#TCGA_enriched_Disease2<-TCGA_enriched_Disease[,-1]

TCGA_geneXenriched_Disease2<-TCGA_geneXenriched_Disease[order(match(TCGA_geneXenriched_Disease$Gene, TCGA_geneXenriched_GOBP$Gene)),]

write_tsv(TCGA_geneXenriched_Disease2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_GeneXenriched_Disease.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de Disease en funcion de TCGA_geneXenriched_Disease$Gene

cormatrix_TCGA_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_Disease_corMatrix.txt", col_names=T)

cormatrix_TCGA_Disease2<-cormatrix_TCGA_Disease[order(match(cormatrix_TCGA_Disease$PE1, TCGA_geneXenriched_Disease2$Gene)),]

write.table(cormatrix_TCGA_Disease2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_Disease_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_TCGA_Disease3<-cormatrix_TCGA_Disease2[,-1]
#sum(colnames(TCGA_enriched_GOBP2) == colnames(cormatrix_TCGA_Disease3))
#1297
write.table(colnames(cormatrix_TCGA_Disease3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_Disease_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#### ordenar TCGA MSIGDB

TCGA_enriched_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_MSigDB/TCGA_enriched_MSigDB.txt", col_names=T)

TCGA_geneXenriched_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_MSigDB/TCGA_GeneXenriched_MSigDB.txt", col_names=T)

#TCGA_enriched_MSigDB2<-TCGA_enriched_MSigDB[,-1]

TCGA_geneXenriched_MSigDB2<-TCGA_geneXenriched_MSigDB[order(match(TCGA_geneXenriched_MSigDB$Gene, TCGA_geneXenriched_GOBP$Gene)),]

write_tsv(TCGA_geneXenriched_MSigDB2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_MSigDB/TCGA_GeneXenriched_MSigDB.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de MSigDB en funcion de TCGA_geneXenriched_MSigDB$Gene

cormatrix_TCGA_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_MSigDB/TCGA_MSigDB_corMatrix.txt")

cormatrix_TCGA_MSigDB2<-cormatrix_TCGA_MSigDB[order(match(cormatrix_TCGA_MSigDB$PE1, TCGA_geneXenriched_MSigDB2$Gene)),]

write.table(cormatrix_TCGA_MSigDB2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_MSigDB/TCGA_MSigDB_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_TCGA_MSigDB3<-cormatrix_TCGA_MSigDB2[,-1]
#sum(colnames(TCGA_enriched_GOBP2) == colnames(cormatrix_TCGA_MSigDB3))
#1297
write.table(colnames(cormatrix_TCGA_MSigDB3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_MSigDB/TCGA_MSigDB_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")


##############################################################################################################

#CCLE CCLE CCLE CCLE


##################################################################################################################33


CCLE_enriched_GOBP<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_enriched_GOBP.txt", col_names=T)

CCLE_geneXenriched_GOBP<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_GeneXenriched_GOBP.txt", col_names=T)
#sum(colnames(CCLE_geneXenriched_GOBP2) == CCLE_enriched_GOBP$GO_id)
#[1] 11854
#> dim(CCLE_geneXenriched_GOBP2)
#[1] 15296 11854
#> dim(CCLE_enriched_GOBP2)
#[1] 11854  1297

cormatrix<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_GOBP_corMatrix.txt")
#> sum(cormatrix$PE1 == CCLE_geneXenriched_GOBP$Gene)

##### GOCC ###############################################3
### ahora hay que ordenar GOCC en funcion de los names que ha generado  GOBP

CCLE_enriched_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_enriched_GOCC.txt", col_names=T)

CCLE_geneXenriched_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_GeneXenriched_GOCC.txt", col_names=T)

#CCLE_enriched_GOCC2<-CCLE_enriched_GOCC[,-1]
# dim(CCLE_enriched_GOCC2)
#[1] 1339 1297
# head(CCLE_enriched_GOCC2[1:2,1:2])
# A tibble: 2 x 2
#  ENSG00000000460 ENSG00000003249
#            <dbl>           <dbl>
#1               0               0
#2               0               0
# sum(colnames(CCLE_enriched_GOCC2) == colnames(CCLE_enriched_GOBP2))
#[1] 1297
# sum(CCLE_geneXenriched_GOBP$Gene == CCLE_geneXenriched_GOCC$Gene)
### 172 (hay que reordenar)
CCLE_geneXenriched_GOCC2<-CCLE_geneXenriched_GOCC[order(match(CCLE_geneXenriched_GOCC$Gene, CCLE_geneXenriched_GOBP$Gene)),]

write_tsv(CCLE_geneXenriched_GOCC2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_GeneXenriched_GOCC.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de GOCC en funcion de CCLE_geneXenriched_GOCC$Gene

cormatrix_CCLE_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_GOCC_corMatrix.txt", col_names=T)
#sum(cormatrix_CCLE_GOCC$PE1 == CCLE_geneXenriched_GOCC2$Gene)
#172
cormatrix_CCLE_GOCC2<-cormatrix_CCLE_GOCC[order(match(cormatrix_CCLE_GOCC$PE1, CCLE_geneXenriched_GOCC2$Gene)),]
#sum(cormatrix_CCLE_GOCC2$PE1 == CCLE_geneXenriched_GOCC2$Gene)
write.table(cormatrix_CCLE_GOCC2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_GOCC_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(cormatrix_CCLE_GOCC2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")


#comprobar los uPE (colnames)
cormatrix_CCLE_GOCC3<-cormatrix_CCLE_GOCC2[,-1]
#sum(colnames(CCLE_enriched_GOBP2 == colnames(cormatrix_CCLE_GOCC3)))
#1297
write.table(colnames(cormatrix_CCLE_GOCC3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_GOCC_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(colnames(cormatrix_CCLE_GOCC3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

##### ORDENAR CCLE_C y CCLE_corFDR para que tambien tengan el mismo orden

CCLE_c<-get(load("~/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_c.Rdata"))
CCLE_c2<-CCLE_c[order(match(rownames(CCLE_c), cormatrix_CCLE_GOCC2$PE1)),]
#sum(colnames(CCLE_enriched_GOCC)== colnames(CCLE_c2))
save(CCLE_c2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_c.Rdata")

CCLE_cor<-get(load("~/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_cor.Rdata"))
CCLE_cor2<-CCLE_cor[order(match(rownames(CCLE_cor), cormatrix_CCLE_GOCC2$PE1)),]
save(CCLE_cor2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_cor.Rdata")

CCLE_corP<-get(load("~/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_corP.Rdata"))
CCLE_corP2<-CCLE_corP[order(match(rownames(CCLE_corP), cormatrix_CCLE_GOCC2$PE1)),]
save(CCLE_corP2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_corP.Rdata")

CCLE_corFDR<-get(load("~/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_corFDR.Rdata"))
CCLE_corFDR2<-CCLE_corFDR[order(match(rownames(CCLE_corFDR), cormatrix_CCLE_GOCC2$PE1)),]
save(CCLE_corFDR2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/CCLE_corFDR.Rdata")

####### ORDENAR GOMF

CCLE_enriched_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_enriched_GOMF.txt", col_names=T)

CCLE_geneXenriched_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_GeneXenriched_GOMF.txt", col_names=T)

#CCLE_enriched_GOMF2<-CCLE_enriched_GOMF[,-1]

CCLE_geneXenriched_GOMF2<-CCLE_geneXenriched_GOMF[order(match(CCLE_geneXenriched_GOMF$Gene, CCLE_geneXenriched_GOBP$Gene)),]

write_tsv(CCLE_geneXenriched_GOMF2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_GeneXenriched_GOMF.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de GOMF en funcion de CCLE_geneXenriched_GOMF$Gene

cormatrix_CCLE_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_GOMF_corMatrix.txt")

cormatrix_CCLE_GOMF2<-cormatrix_CCLE_GOMF[order(match(cormatrix_CCLE_GOMF$PE1, CCLE_geneXenriched_GOMF2$Gene)),]

write.table(cormatrix_CCLE_GOMF2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_GOMF_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_CCLE_GOMF3<-cormatrix_CCLE_GOMF2[,-1]
#sum(colnames(CCLE_enriched_GOBP2) == colnames(cormatrix_CCLE_GOMF3))
#1297
write.table(colnames(cormatrix_CCLE_GOMF3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_GOMF_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

####### ORDENAR Disease

CCLE_enriched_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_enriched_Disease.txt", col_names=T)

CCLE_geneXenriched_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_GeneXenriched_Disease.txt", col_names=T)

#CCLE_enriched_Disease2<-CCLE_enriched_Disease[,-1]

CCLE_geneXenriched_Disease2<-CCLE_geneXenriched_Disease[order(match(CCLE_geneXenriched_Disease$Gene, CCLE_geneXenriched_GOBP$Gene)),]

write_tsv(CCLE_geneXenriched_Disease2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_GeneXenriched_Disease.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de Disease en funcion de CCLE_geneXenriched_Disease$Gene

cormatrix_CCLE_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_Disease_corMatrix.txt", col_names=T)

cormatrix_CCLE_Disease2<-cormatrix_CCLE_Disease[order(match(cormatrix_CCLE_Disease$PE1, CCLE_geneXenriched_Disease2$Gene)),]

write.table(cormatrix_CCLE_Disease2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_Disease_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_CCLE_Disease3<-cormatrix_CCLE_Disease2[,-1]
#sum(colnames(CCLE_enriched_GOBP2) == colnames(cormatrix_CCLE_Disease3))
#1297
write.table(colnames(cormatrix_CCLE_Disease3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_Disease_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#### ordenar CCLE MSIGDB

CCLE_enriched_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_MSigDB/CCLE_enriched_MSigDB.txt", col_names=T)

CCLE_geneXenriched_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_MSigDB/CCLE_GeneXenriched_MSigDB.txt", col_names=T)

#CCLE_enriched_MSigDB2<-CCLE_enriched_MSigDB[,-1]

CCLE_geneXenriched_MSigDB2<-CCLE_geneXenriched_MSigDB[order(match(CCLE_geneXenriched_MSigDB$Gene, CCLE_geneXenriched_GOBP$Gene)),]

write_tsv(CCLE_geneXenriched_MSigDB2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_MSigDB/CCLE_GeneXenriched_MSigDB.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de MSigDB en funcion de CCLE_geneXenriched_MSigDB$Gene

cormatrix_CCLE_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_MSigDB/CCLE_MSigDB_corMatrix.txt")

cormatrix_CCLE_MSigDB2<-cormatrix_CCLE_MSigDB[order(match(cormatrix_CCLE_MSigDB$PE1, CCLE_geneXenriched_MSigDB2$Gene)),]

write.table(cormatrix_CCLE_MSigDB2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_MSigDB/CCLE_MSigDB_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_CCLE_MSigDB3<-cormatrix_CCLE_MSigDB2[,-1]
#sum(colnames(CCLE_enriched_GOBP2) == colnames(cormatrix_CCLE_MSigDB3))
#1297
write.table(colnames(cormatrix_CCLE_MSigDB3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_MSigDB/CCLE_MSigDB_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")


#####################################################

GTEX GTEX GTEX

####################################################3

GTEX_enriched_GOBP<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_enriched_GOBP.txt", col_names=T)

GTEX_geneXenriched_GOBP<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_GeneXenriched_GOBP.txt", col_names=T)
#sum(colnames(GTEX_geneXenriched_GOBP2) == GTEX_enriched_GOBP$GO_id)
#[1] 11854
#> dim(GTEX_geneXenriched_GOBP2)
#[1] 15296 11854
#> dim(GTEX_enriched_GOBP2)
#[1] 11854  1297

cormatrix<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_GOBP_corMatrix.txt")
#> sum(cormatrix$PE1 == GTEX_geneXenriched_GOBP$Gene)

##### GOCC ###############################################3
### ahora hay que ordenar GOCC en funcion de los names que ha generado  GOBP

GTEX_enriched_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_enriched_GOCC.txt", col_names=T)

GTEX_geneXenriched_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_GeneXenriched_GOCC.txt", col_names=T)

#GTEX_enriched_GOCC2<-GTEX_enriched_GOCC[,-1]
# dim(GTEX_enriched_GOCC2)
#[1] 1339 1297
# head(GTEX_enriched_GOCC2[1:2,1:2])
# A tibble: 2 x 2
#  ENSG00000000460 ENSG00000003249
#            <dbl>           <dbl>
#1               0               0
#2               0               0
# sum(colnames(GTEX_enriched_GOCC2) == colnames(GTEX_enriched_GOBP2))
#[1] 1297
# sum(GTEX_geneXenriched_GOBP$Gene == GTEX_geneXenriched_GOCC$Gene)
### 172 (hay que reordenar)
GTEX_geneXenriched_GOCC2<-GTEX_geneXenriched_GOCC[order(match(GTEX_geneXenriched_GOCC$Gene, GTEX_geneXenriched_GOBP$Gene)),]

write_tsv(GTEX_geneXenriched_GOCC2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_GeneXenriched_GOCC.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de GOCC en funcion de GTEX_geneXenriched_GOCC$Gene

cormatrix_GTEX_GOCC<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_GOCC_corMatrix.txt", col_names=T)
#sum(cormatrix_GTEX_GOCC$PE1 == GTEX_geneXenriched_GOCC2$Gene)
#172
cormatrix_GTEX_GOCC2<-cormatrix_GTEX_GOCC[order(match(cormatrix_GTEX_GOCC$PE1, GTEX_geneXenriched_GOCC2$Gene)),]
#sum(cormatrix_GTEX_GOCC2$PE1 == GTEX_geneXenriched_GOCC2$Gene)
write.table(cormatrix_GTEX_GOCC2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_GOCC_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(cormatrix_GTEX_GOCC2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")


#comprobar los uPE (colnames)
cormatrix_GTEX_GOCC3<-cormatrix_GTEX_GOCC2[,-1]
#sum(colnames(GTEX_enriched_GOBP2 == colnames(cormatrix_GTEX_GOCC3)))
#1297
write.table(colnames(cormatrix_GTEX_GOCC3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_GOCC_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(colnames(cormatrix_GTEX_GOCC3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

##### ORDENAR GTEX_C y GTEX_corFDR para que tambien tengan el mismo orden

GTEX_c<-get(load("~/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_c.Rdata"))
GTEX_c2<-GTEX_c[order(match(rownames(GTEX_c), cormatrix_GTEX_GOCC2$PE1)),]
sum(colnames(GTEX_enriched_GOBP2)== colnames(GTEX_c2))
save(GTEX_c2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_c.Rdata")

GTEX_cor<-get(load("~/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_cor.Rdata"))
GTEX_cor2<-GTEX_cor[order(match(rownames(GTEX_cor), cormatrix_GTEX_GOCC2$PE1)),]
save(GTEX_cor2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_cor.Rdata")

GTEX_corP<-get(load("~/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_corP.Rdata"))
GTEX_corP2<-GTEX_corP[order(match(rownames(GTEX_corP), cormatrix_GTEX_GOCC2$PE1)),]
save(GTEX_corP2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_corP.Rdata")

GTEX_corFDR<-get(load("~/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_corFDR.Rdata"))
GTEX_corFDR2<-GTEX_corFDR[order(match(rownames(GTEX_corFDR), cormatrix_GTEX_GOCC2$PE1)),]
save(GTEX_corFDR2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/GTEX_corFDR.Rdata")

####### ORDENAR GOMF

#GTEX_enriched_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_enriched_GOMF.txt", col_names=T)

GTEX_geneXenriched_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_GeneXenriched_GOMF.txt", col_names=T)

#GTEX_enriched_GOMF2<-GTEX_enriched_GOMF[,-1]

GTEX_geneXenriched_GOMF2<-GTEX_geneXenriched_GOMF[order(match(GTEX_geneXenriched_GOMF$Gene, GTEX_geneXenriched_GOBP$Gene)),]

write_tsv(GTEX_geneXenriched_GOMF2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_GeneXenriched_GOMF.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de GOMF en funcion de GTEX_geneXenriched_GOMF$Gene

cormatrix_GTEX_GOMF<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_GOMF_corMatrix.txt")

cormatrix_GTEX_GOMF2<-cormatrix_GTEX_GOMF[order(match(cormatrix_GTEX_GOMF$PE1, GTEX_geneXenriched_GOMF2$Gene)),]

write.table(cormatrix_GTEX_GOMF2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_GOMF_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_GTEX_GOMF3<-cormatrix_GTEX_GOMF2[,-1]
#sum(colnames(GTEX_enriched_GOBP2) == colnames(cormatrix_GTEX_GOMF3))
#1297
write.table(colnames(cormatrix_GTEX_GOMF3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_GOMF_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

####### ORDENAR Disease

GTEX_enriched_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_enriched_Disease.txt", col_names=T)

GTEX_geneXenriched_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_GeneXenriched_Disease.txt", col_names=T)

#GTEX_enriched_Disease2<-GTEX_enriched_Disease[,-1]

GTEX_geneXenriched_Disease2<-GTEX_geneXenriched_Disease[order(match(GTEX_geneXenriched_Disease$Gene, GTEX_geneXenriched_GOBP$Gene)),]

write_tsv(GTEX_geneXenriched_Disease2, path="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_GeneXenriched_Disease.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de Disease en funcion de GTEX_geneXenriched_Disease$Gene

cormatrix_GTEX_Disease<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_Disease_corMatrix.txt", col_names=T)

cormatrix_GTEX_Disease2<-cormatrix_GTEX_Disease[order(match(cormatrix_GTEX_Disease$PE1, GTEX_geneXenriched_Disease2$Gene)),]

write.table(cormatrix_GTEX_Disease2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_Disease_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_GTEX_Disease3<-cormatrix_GTEX_Disease2[,-1]
#sum(colnames(GTEX_enriched_GOBP2) == colnames(cormatrix_GTEX_Disease3))
#1297
write.table(colnames(cormatrix_GTEX_Disease3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_Disease_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#### ordenar GTEX MSIGDB

GTEX_enriched_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_MSigDB/GTEX_enriched_MSigDB.txt", col_names=T)

GTEX_geneXenriched_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_MSigDB/GTEX_GeneXenriched_MSigDB.txt", col_names=T)

#GTEX_enriched_MSigDB2<-GTEX_enriched_MSigDB[,-1]

GTEX_geneXenriched_MSigDB2<-GTEX_geneXenriched_MSigDB[order(match(GTEX_geneXenriched_MSigDB$Gene, GTEX_geneXenriched_GOBP$Gene)),]

write_tsv(GTEX_geneXenriched_MSigDB2, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_MSigDB/GTEX_GeneXenriched_MSigDB.txt", col_names=T, quote=F)
 
## ahora tengo que ordenar la cormatrix$PE1 de MSigDB en funcion de GTEX_geneXenriched_MSigDB$Gene

cormatrix_GTEX_MSigDB<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_MSigDB/GTEX_MSigDB_corMatrix.txt")

cormatrix_GTEX_MSigDB2<-cormatrix_GTEX_MSigDB[order(match(cormatrix_GTEX_MSigDB$PE1, GTEX_geneXenriched_MSigDB2$Gene)),]

write.table(cormatrix_GTEX_MSigDB2$PE1, file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_MSigDB/GTEX_MSigDB_PE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

#comprobar los uPE (colnames)
cormatrix_GTEX_MSigDB3<-cormatrix_GTEX_MSigDB2[,-1]
#sum(colnames(GTEX_enriched_GOBP2) == colnames(cormatrix_GTEX_MSigDB3))
#1297
write.table(colnames(cormatrix_GTEX_MSigDB3), file="/home/nostromo/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_MSigDB/GTEX_MSigDB_UPE_names.txt", col.names=F, row.names=F, quote=F, sep="\t")

