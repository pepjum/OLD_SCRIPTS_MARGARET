########################################################33

##### TCGA

###########################################################


## GOBP
library(tidyverse)

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOBP_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOBP_unselected_GOS<-paste(GOS_unselected)

TCGA_enriched_GOBP<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_enriched_GOBP.txt", col_names=T)

TCGA_enriched_GOBP_filtered<-TCGA_enriched_GOBP[-which(TCGA_enriched_GOBP$GO_id %in% GOBP_unselected_GOS),]
write_tsvTCGA_enriched_GOBP_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_enriched_GOBP_filtered.txt", col_names=T, quote=T)

TCGA_geneXenriched_GOBP<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_GeneXenriched_GOBP.txt", col_names=T)
TCGA_geneXenriched_GOBP_filtered<-TCGA_geneXenriched_GOBP[,-which(colnames(TCGA_geneXenriched_GOBP) %in% GOBP_unselected_GOS)]
write_tsv(TCGA_geneXenriched_GOBP_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_GeneXenriched_GOBP_filtered.txt", col_names=T, quote=F)

write.table(TCGA_enriched_GOBP_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_GOBP_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

TCGA_enriched_GOBP_PVAL<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_enriched_GOBP_PVAL.txt", col_names=T)

TCGA_enriched_GOBP_PVAL_filtered<-TCGA_enriched_GOBP_PVAL[-which(TCGA_enriched_GOBP_PVAL$GO_id %in% GOBP_unselected_GOS),]
write.table(TCGA_enriched_GOBP_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOBP/TCGA_enriched_GOBP_PVAL_filtered.txt", col_names=T, quote=T)


#### GOCC

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOCC_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOCC_unselected_GOS<-paste(GOS_unselected)

TCGA_enriched_GOCC<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_enriched_GOCC.txt", col_names=T)

TCGA_enriched_GOCC_filtered<-TCGA_enriched_GOCC[-which(TCGA_enriched_GOCC$GO_id %in% GOCC_unselected_GOS),]
write_tsv(TCGA_enriched_GOCC_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_enriched_GOCC_filtered.txt", col_names=T, quote=F)

TCGA_geneXenriched_GOCC<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_GeneXenriched_GOCC.txt", col_names=T)
TCGA_geneXenriched_GOCC_filtered<-TCGA_geneXenriched_GOCC[,-which(colnames(TCGA_geneXenriched_GOCC) %in% GOCC_unselected_GOS)]
write_tsv(TCGA_geneXenriched_GOCC_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_GeneXenriched_GOCC_filtered.txt", col_names=T, quote=F)

write.table(TCGA_enriched_GOCC_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_GOCC_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

TCGA_enriched_GOCC_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_enriched_GOCC_PVAL.txt", header=T, check.names=F)

TCGA_enriched_GOCC_PVAL_filtered<-TCGA_enriched_GOCC_PVAL[-which(TCGA_enriched_GOCC_PVAL$GO_id %in% GOCC_unselected_GOS),]
write.table(TCGA_enriched_GOCC_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOCC/TCGA_enriched_GOCC_PVAL_filtered.txt", col.names=T, quote=T, sep="\t", row.names=F)



#### GOMF

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOMF_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOMF_unselected_GOS<-paste(GOS_unselected)

TCGA_enriched_GOMF<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_enriched_GOMF.txt", col_names=T)

TCGA_enriched_GOMF_filtered<-TCGA_enriched_GOMF[-which(TCGA_enriched_GOMF$GO_id %in% GOMF_unselected_GOS),]
write_tsv(TCGA_enriched_GOMF_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_enriched_GOMF_filtered.txt", col_names=T, quote=F)

TCGA_geneXenriched_GOMF<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_GeneXenriched_GOMF.txt", col_names=T)
TCGA_geneXenriched_GOMF_filtered<-TCGA_geneXenriched_GOMF[,-which(colnames(TCGA_geneXenriched_GOMF) %in% GOMF_unselected_GOS)]
write_tsv(TCGA_geneXenriched_GOMF_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_GeneXenriched_GOMF_filtered.txt", col_names=T, quote=F)

write.table(TCGA_enriched_GOMF_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_GOMF_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

TCGA_enriched_GOMF_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_enriched_GOMF_PVAL.txt", header=T, check.names=F)

TCGA_enriched_GOMF_PVAL_filtered<-TCGA_enriched_GOMF_PVAL[-which(TCGA_enriched_GOMF_PVAL$GO_id %in% GOMF_unselected_GOS),]
write.table(TCGA_enriched_GOMF_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_GOMF/TCGA_enriched_GOMF_PVAL_filtered.txt", col.names=T, quote=T, sep="\t", row.names=F)


#### Disease

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/Disease_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

Disease_unselected_GOS<-paste(GOS_unselected)

TCGA_enriched_Disease<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_enriched_Disease.txt", col_names=T)

TCGA_enriched_Disease_filtered<-TCGA_enriched_Disease[-which(TCGA_enriched_Disease$GO_id %in% Disease_unselected_GOS),]
write_tsv(TCGA_enriched_Disease_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_enriched_Disease_filtered.txt", col_names=T, quote=F)

TCGA_enriched_Disease_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_enriched_Disease_PVAL.txt", header=T, check.names=F)

TCGA_enriched_Disease_PVAL_filtered<-TCGA_enriched_Disease_PVAL[-which(TCGA_enriched_Disease_PVAL$GO_id %in% Disease_unselected_GOS),]
write.table(TCGA_enriched_Disease_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_enriched_Disease_PVAL_filtered.txt", col.names=T, quote=T, sep="\t", row.names=F)



TCGA_geneXenriched_Disease<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_GeneXenriched_Disease.txt", col_names=T)
TCGA_geneXenriched_Disease_filtered<-TCGA_geneXenriched_Disease[,-which(colnames(TCGA_geneXenriched_Disease) %in% Disease_unselected_GOS)]
write_tsv(TCGA_geneXenriched_Disease_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_GeneXenriched_Disease_filtered.txt", col_names=T, quote=F)

write.table(TCGA_enriched_Disease_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/TCGA/OUTPUT_TCGA_Disease/TCGA_Disease_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)


########################################################33

##### GTEX

###########################################################


## GOBP
library(tidyverse)

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOBP_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOBP_unselected_GOS<-paste(GOS_unselected)

GTEX_enriched_GOBP<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_enriched_GOBP.txt", col_names=T)

GTEX_enriched_GOBP_filtered<-GTEX_enriched_GOBP[-which(GTEX_enriched_GOBP$GO_id %in% GOBP_unselected_GOS),]
write_tsv(GTEX_enriched_GOBP_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_enriched_GOBP_filtered.txt", col_names=T, quote=F)

GTEX_geneXenriched_GOBP<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_GeneXenriched_GOBP.txt", col_names=T)
GTEX_geneXenriched_GOBP_filtered<-GTEX_geneXenriched_GOBP[,-which(colnames(GTEX_geneXenriched_GOBP) %in% GOBP_unselected_GOS)]
write_tsv(GTEX_geneXenriched_GOBP_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_GeneXenriched_GOBP_filtered.txt", col_names=T, quote=F)

write.table(GTEX_enriched_GOBP_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_GOBP_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

GTEX_enriched_GOBP_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_enriched_PVAL_GOBP.txt", header=T, check.names=F)

GTEX_enriched_GOBP_PVAL_filtered<-GTEX_enriched_GOBP_PVAL[-which(GTEX_enriched_GOBP_PVAL$GO_id %in% GOBP_unselected_GOS),]
write.table(GTEX_enriched_GOBP_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOBP/GTEX_enriched_GOBP_PVAL_filtered.txt", col.names=T, quote=T, sep="\t", row.names=F)


#### GOCC

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOCC_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOCC_unselected_GOS<-paste(GOS_unselected)

GTEX_enriched_GOCC<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_enriched_GOCC.txt", col_names=T)

GTEX_enriched_GOCC_filtered<-GTEX_enriched_GOCC[-which(GTEX_enriched_GOCC$GO_id %in% GOCC_unselected_GOS),]
write_tsv(GTEX_enriched_GOCC_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_enriched_GOCC_filtered.txt", col_names=T, quote=F)

GTEX_geneXenriched_GOCC<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_GeneXenriched_GOCC.txt", col_names=T)
GTEX_geneXenriched_GOCC_filtered<-GTEX_geneXenriched_GOCC[,-which(colnames(GTEX_geneXenriched_GOCC) %in% GOCC_unselected_GOS)]
write_tsv(GTEX_geneXenriched_GOCC_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_GeneXenriched_GOCC_filtered.txt", col_names=T, quote=F)

write.table(GTEX_enriched_GOCC_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_GOCC_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

GTEX_enriched_GOCC_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_enriched_GOCC_PVAL.txt", header=T)

GTEX_enriched_GOCC_PVAL_filtered<-GTEX_enriched_GOCC_PVAL[-which(GTEX_enriched_GOCC_PVAL$GO_id %in% GOCC_unselected_GOS),]
write.table(GTEX_enriched_GOCC_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOCC/GTEX_enriched_GOCC_PVAL_filtered.txt", col.names=T, quote=T, row.names=F, sep="\t")



#### GOMF

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOMF_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOMF_unselected_GOS<-paste(GOS_unselected)

GTEX_enriched_GOMF<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_enriched_GOMF.txt", col_names=T)

GTEX_enriched_GOMF_filtered<-GTEX_enriched_GOMF[-which(GTEX_enriched_GOMF$GO_id %in% GOMF_unselected_GOS),]
write_tsv(GTEX_enriched_GOMF_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_enriched_GOMF_filtered.txt", col_names=T, quote=F)

GTEX_geneXenriched_GOMF<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_GeneXenriched_GOMF.txt", col_names=T)
GTEX_geneXenriched_GOMF_filtered<-GTEX_geneXenriched_GOMF[,-which(colnames(GTEX_geneXenriched_GOMF) %in% GOMF_unselected_GOS)]
write_tsv(GTEX_geneXenriched_GOMF_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_GeneXenriched_GOMF_filtered.txt", col_names=T, quote=F)

write.table(GTEX_enriched_GOMF_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_GOMF_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

GTEX_enriched_GOMF_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_enriched_GOMF_PVAL.txt", header=T, check.names=F)

GTEX_enriched_GOMF_PVAL_filtered<-GTEX_enriched_GOMF_PVAL[-which(GTEX_enriched_GOMF_PVAL$GO_id %in% GOMF_unselected_GOS),]
write.table(GTEX_enriched_GOMF_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_GOMF/GTEX_enriched_GOMF_PVAL_filtered.txt", col.names=T, quote=T, sep="\t", row.names=F)


#### Disease

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/Disease_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

Disease_unselected_GOS<-paste(GOS_unselected)

GTEX_enriched_Disease<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_enriched_Disease.txt", col_names=T)

GTEX_enriched_Disease_filtered<-GTEX_enriched_Disease[-which(GTEX_enriched_Disease$GO_id %in% Disease_unselected_GOS),]
write_tsv(GTEX_enriched_Disease_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_enriched_Disease_filtered.txt", col_names=T, quote=F)

GTEX_geneXenriched_Disease<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_GeneXenriched_Disease.txt", col_names=T)
GTEX_geneXenriched_Disease_filtered<-GTEX_geneXenriched_Disease[,-which(colnames(GTEX_geneXenriched_Disease) %in% Disease_unselected_GOS)]
write_tsv(GTEX_geneXenriched_Disease_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_GeneXenriched_Disease_filtered.txt", col_names=T, quote=F)

write.table(GTEX_enriched_Disease_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_Disease_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

GTEX_enriched_Disease_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_enriched_Disease_PVAL.txt", header=T, check.names=F)

GTEX_enriched_Disease_PVAL_filtered<-GTEX_enriched_Disease_PVAL[-which(GTEX_enriched_Disease_PVAL$GO_id %in% Disease_unselected_GOS),]
write.table(GTEX_enriched_Disease_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/GTEX/OUTPUT_GTEX_Disease/GTEX_enriched_Disease_PVAL_filtered.txt", col.names=T,row.names=F, quote=T, sep="\t")




########################################################33

##### CCLE

###########################################################


## GOBP
library(tidyverse)

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOBP_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOBP_unselected_GOS<-paste(GOS_unselected)

CCLE_enriched_GOBP<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_enriched_GOBP.txt", col_names=T)

CCLE_enriched_GOBP_filtered<-CCLE_enriched_GOBP[-which(CCLE_enriched_GOBP$GO_id %in% GOBP_unselected_GOS),]
write_tsv(CCLE_enriched_GOBP_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_enriched_GOBP_filtered.txt", col_names=T, quote=F)

CCLE_geneXenriched_GOBP<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_GeneXenriched_GOBP.txt", col_names=T)
CCLE_geneXenriched_GOBP_filtered<-CCLE_geneXenriched_GOBP[,-which(colnames(CCLE_geneXenriched_GOBP) %in% GOBP_unselected_GOS)]
write_tsv(CCLE_geneXenriched_GOBP_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_GeneXenriched_GOBP_filtered.txt", col_names=T, quote=F)

write.table(CCLE_enriched_GOBP_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_GOBP_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

CCLE_enriched_GOBP_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_enriched_GOBP_PVAL.txt", header=T, check.names=F)

CCLE_enriched_GOBP_PVAL_filtered<-CCLE_enriched_GOBP_PVAL[-which(CCLE_enriched_GOBP_PVAL$GO_id %in% GOBP_unselected_GOS),]
write.table(CCLE_enriched_GOBP_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOBP/CCLE_enriched_GOBP_PVAL_filtered.txt", col.names=T, quote=T, row.names=F, sep="\t")


#### GOCC

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOCC_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOCC_unselected_GOS<-paste(GOS_unselected)

CCLE_enriched_GOCC<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_enriched_GOCC.txt", col_names=T)

CCLE_enriched_GOCC_filtered<-CCLE_enriched_GOCC[-which(CCLE_enriched_GOCC$GO_id %in% GOCC_unselected_GOS),]
write_tsv(CCLE_enriched_GOCC_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_enriched_GOCC_filtered.txt", col_names=T, quote=F)

CCLE_geneXenriched_GOCC<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_GeneXenriched_GOCC.txt", col_names=T)
CCLE_geneXenriched_GOCC_filtered<-CCLE_geneXenriched_GOCC[,-which(colnames(CCLE_geneXenriched_GOCC) %in% GOCC_unselected_GOS)]
write_tsv(CCLE_geneXenriched_GOCC_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_GeneXenriched_GOCC_filtered.txt", col_names=T, quote=F)

write.table(CCLE_enriched_GOCC_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_GOCC_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

CCLE_enriched_GOCC_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_enriched_GOCC_PVAL.txt", header=T, check.names=F)

CCLE_enriched_GOCC_PVAL_filtered<-CCLE_enriched_GOCC_PVAL[-which(CCLE_enriched_GOCC_PVAL$GO_id %in% GOCC_unselected_GOS),]
write.table(CCLE_enriched_GOCC_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOCC/CCLE_enriched_GOCC_PVAL_filtered.txt", col.names=T, quote=T, sep="\t", row.names=F)



#### GOMF

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/GOMF_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

GOMF_unselected_GOS<-paste(GOS_unselected)

CCLE_enriched_GOMF<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_enriched_GOMF.txt", col_names=T)

CCLE_enriched_GOMF_filtered<-CCLE_enriched_GOMF[-which(CCLE_enriched_GOMF$GO_id %in% GOMF_unselected_GOS),]
write_tsv(CCLE_enriched_GOMF_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_enriched_GOMF_filtered.txt", col_names=T, quote=F)

CCLE_geneXenriched_GOMF<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_GeneXenriched_GOMF.txt", col_names=T)
CCLE_geneXenriched_GOMF_filtered<-CCLE_geneXenriched_GOMF[,-which(colnames(CCLE_geneXenriched_GOMF) %in% GOMF_unselected_GOS)]
write_tsv(CCLE_geneXenriched_GOMF_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_GeneXenriched_GOMF_filtered.txt", col_names=T, quote=F)

write.table(CCLE_enriched_GOMF_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_GOMF_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

CCLE_enriched_GOMF_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_enriched_GOMF_PVAL.txt", header=T, check.names=F)

CCLE_enriched_GOMF_PVAL_filtered<-CCLE_enriched_GOMF_PVAL[-which(CCLE_enriched_GOMF_PVAL$GO_id %in% GOMF_unselected_GOS),]
write.table(CCLE_enriched_GOMF_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_GOMF/CCLE_enriched_GOMF_PVAL_filtered.txt", col.names=T, row.names=F, quote=T, sep="\t")


#### Disease

GOS_unselected<-read_tsv("/home/nostromo/data/pepe/43_UPEFINDER_AGO20/BlackList/Disease_unselected.txt", col_names=F)
GOS_unselected<-paste(GOS_unselected$X1)
GOS_unselected<-unlist(lapply(strsplit(paste(GOS_unselected)," "),"[",2))

Disease_unselected_GOS<-paste(GOS_unselected)

CCLE_enriched_Disease<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_enriched_Disease.txt", col_names=T)

CCLE_enriched_Disease_filtered<-CCLE_enriched_Disease[-which(CCLE_enriched_Disease$GO_id %in% Disease_unselected_GOS),]
write_tsv(CCLE_enriched_Disease_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_enriched_Disease_filtered.txt", col_names=T, quote=F)

CCLE_geneXenriched_Disease<-read_tsv("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_GeneXenriched_Disease.txt", col_names=T)
CCLE_geneXenriched_Disease_filtered<-CCLE_geneXenriched_Disease[,-which(colnames(CCLE_geneXenriched_Disease) %in% Disease_unselected_GOS)]
write_tsv(CCLE_geneXenriched_Disease_filtered, path="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_GeneXenriched_Disease_filtered.txt", col_names=T, quote=F)

write.table(CCLE_enriched_Disease_filtered$GO_id, file="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_Disease_names_filtered.txt", row.names=F, col.names=F, sep="\t", quote=F)

CCLE_enriched_Disease_PVAL<-read.table("~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_enriched_Disease_PVAL.txt", header=T, check.names=F)

CCLE_enriched_Disease_PVAL_filtered<-CCLE_enriched_Disease_PVAL[-which(CCLE_enriched_Disease_PVAL$GO_id %in% Disease_unselected_GOS),]
write.table(CCLE_enriched_Disease_PVAL_filtered, file="~/data/pepe/43_UPEFINDER_AGO20/CCLE/OUTPUT_CCLE_Disease/CCLE_enriched_Disease_PVAL_filtered.txt", col.names=T, row.names=F, quote=T, sep="\t")

