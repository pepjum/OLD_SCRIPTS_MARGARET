setwd("~/dato-activo/03_Analysis/eguruce/20_UPs_Abr19/Experiments")
## Sevastopol
setwd("~/data/eguruce/20_UPs_Abr19")


# save.image("UPsVal_Feb20.RData")
# load("UPsVal_Feb20.RData")

# save.image("UPsVal_paper_Ago20.RData")
# load("UPsVal_paper_Ago20.RData")


## All TERMS x PageRanks
pageRank_GO_CCLE <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_GO_CCLE_Aug20.tsv",header=F, sep="\t")
pageRank_Dis_CCLE <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_Mala_CCLE_Aug20.tsv",header=F, sep="\t")
pageRank_MSig_CCLE <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_MSIG_CCLE_Aug20.tsv",header=F, sep="\t")
pageRank_GO_GTEX <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_GO_GTEX_Aug20.tsv",header=F, sep="\t")
pageRank_Dis_GTEX <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_Mala_GTEX_Aug20.tsv",header=F, sep="\t")
pageRank_MSig_GTEX <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_MSIG_GTEX_Aug20.tsv",header=F, sep="\t")
pageRank_GO_TCGA <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_GO_TCGA_Aug20.tsv",header=F, sep="\t")
pageRank_Dis_TCGA <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_Mala_TCGA_Aug20.tsv",header=F, sep="\t")
pageRank_MSig_TCGA <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_MSIG_TCGA_Aug20.tsv",header=F, sep="\t")
GOdb <- read.table("GOannotation_Ago20.txt", header=T, sep="\t")
rownames(GOdb) <- GOdb[,1]

DB_ggplot2 <- rbind(data.frame(pageRank_GO_CCLE[paste(pageRank_GO_CCLE[,1]) %in% rownames(uPExGOBP_CCLE),],DB="GO_BP",Experiment="CCLE"),
  data.frame(pageRank_GO_CCLE[paste(pageRank_GO_CCLE[,1]) %in% rownames(uPExGOMF_CCLE),],DB="GO_MF",Experiment="CCLE"),
 data.frame(pageRank_GO_CCLE[paste(pageRank_GO_CCLE[,1]) %in% rownames(uPExGOCC_CCLE),],DB="GO_CC",Experiment="CCLE"),
 data.frame(pageRank_Dis_CCLE,DB="Disease",Experiment="CCLE"),
 data.frame(pageRank_MSig_CCLE,DB="MSigDB",Experiment="CCLE"),
 data.frame(pageRank_GO_GTEX[paste(pageRank_GO_GTEX[,1]) %in% rownames(uPExGOBP_GTEX),],DB="GO_BP",Experiment="GTEX"),
 data.frame(pageRank_GO_GTEX[paste(pageRank_GO_GTEX[,1]) %in% rownames(uPExGOMF_GTEX),],DB="GO_MF",Experiment="GTEX"),
  data.frame(pageRank_GO_GTEX[paste(pageRank_GO_GTEX[,1]) %in% rownames(uPExGOCC_GTEX),],DB="GO_CC",Experiment="GTEX"),
  data.frame(pageRank_Dis_GTEX,DB="Disease",Experiment="GTEX"),
  data.frame(pageRank_MSig_GTEX,DB="MSigDB",Experiment="GTEX"),
  data.frame(pageRank_GO_TCGA[paste(pageRank_GO_TCGA[,1]) %in% rownames(uPExGOBP_TCGA),],DB="GO_BP",Experiment="TCGA"),
  data.frame(pageRank_GO_TCGA[paste(pageRank_GO_TCGA[,1]) %in% rownames(uPExGOMF_TCGA),],DB="GO_MF",Experiment="TCGA"),
  data.frame(pageRank_GO_TCGA[paste(pageRank_GO_TCGA[,1]) %in% rownames(uPExGOCC_TCGA),],DB="GO_CC",Experiment="TCGA"),
  data.frame(pageRank_Dis_TCGA,DB="Disease",Experiment="TCGA"),
  data.frame(pageRank_MSig_TCGA,DB="MSigDB",Experiment="TCGA"))
DB_ggplot2$log2PageRank <- log2(DB_ggplot2[,2])
DB_ggplot2$log10PageRank <- -10*log10(DB_ggplot2[,2])

p <- ggplot(DB_ggplot2, aes(x=log2PageRank, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")
pdf("PageRanks_DBxExp_histogram.pdf")
p
dev.off()

p <- ggplot(DB_ggplot2, aes(x=log2PageRank)) + geom_density(color="black", fill="white") + facet_grid(Experiment ~ DB) + theme_bw()
pdf("PageRanks_DBxExp_density.pdf")
p
dev.off()
p <- ggplot(DB_ggplot2, aes(x=DB,y=log2PageRank)) + geom_boxplot() + ylab("log2(PRscore)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PageRanks_ProteinxExp_boxplot.pdf")
p
dev.off()
p <- ggplot(DB_ggplot2, aes(x=DB,y=log2PageRank)) + geom_violin() + ylab("log2(PRscore)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PageRanks_ProteinxExp_violin.pdf")
p
dev.off()

# uPEs X PageRanks
pageRank_uPE_CCLE <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_uPE_CCLE_Aug20.tsv",header=F, sep="\t")
pageRank_uPE_GTEX <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_uPE_GTEX_Aug20.tsv",header=F, sep="\t")
pageRank_uPE_TCGA <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_uPE_TCGA_Aug20.tsv",header=F, sep="\t")
# PEs X PageRanks
pageRank_PE_CCLE <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_PE_CCLE_Aug20.tsv",header=F, sep="\t")
pageRank_PE_GTEX <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_PE_GTEX_Aug20.tsv",header=F, sep="\t")
pageRank_PE_TCGA <- read.table("../../gserranos/PageRank/Data/Validation_Eli/Validation_ELI_PE_TCGA_Aug20.tsv",header=F, sep="\t")

Protein_ggplot2 <- rbind(data.frame(pageRank_uPE_CCLE,Protein="uPE",Experiment="CCLE"), data.frame(pageRank_uPE_GTEX,Protein="uPE",Experiment="GTEX"), data.frame(pageRank_uPE_TCGA,Protein="uPE",Experiment="TCGA"), data.frame(pageRank_PE_CCLE,Protein="PE",Experiment="CCLE"), data.frame(pageRank_PE_GTEX,Protein="PE",Experiment="GTEX"), data.frame(pageRank_PE_TCGA,Protein="PE",Experiment="TCGA"))
Protein_ggplot2$log2PageRank <- log2(Protein_ggplot2[,2]) 
p <- ggplot(Protein_ggplot2, aes(x=log2PageRank)) + geom_histogram(color="black", fill="white") + facet_grid(Protein ~ Experiment) + theme_bw()
pdf("PageRanks_ProteinxExp.pdf")
p
dev.off()
p <- ggplot(Protein_ggplot2, aes(x=Experiment,y=log2PageRank)) + geom_boxplot() + ylab("log2(PRscore)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Protein ~ .)
pdf("PageRanks_ProteinxExp_boxplot.pdf")
p
dev.off()


## BEST TERMS 
pageRank_res <- read.table("../../gserranos/PageRank/Data/DataValidation_ELIV2.txt", header=F, sep="\t")
# No me la ha dejado Guille (intento sacarla yo)
pageRank_res <- rbind(cbind(rbind(pageRank_GO_CCLE[paste(pageRank_GO_CCLE[,1]) %in% rownames(uPExGOBP_CCLE),][1:10,],pageRank_GO_CCLE[paste(pageRank_GO_CCLE[,1]) %in% rownames(uPExGOMF_CCLE),][1:10,], pageRank_GO_CCLE[paste(pageRank_GO_CCLE[,1]) %in% rownames(uPExGOCC_CCLE),][1:10,]), Dataset="CCLE"), cbind(rbind(pageRank_GO_TCGA[paste(pageRank_GO_TCGA[,1]) %in% rownames(uPExGOBP_TCGA),][1:10,],pageRank_GO_TCGA[paste(pageRank_GO_TCGA[,1]) %in% rownames(uPExGOMF_TCGA),][1:10,], pageRank_GO_TCGA[paste(pageRank_GO_TCGA[,1]) %in% rownames(uPExGOCC_TCGA),][1:10,]), Dataset="TCGA"), cbind(rbind(pageRank_GO_GTEX[paste(pageRank_GO_GTEX[,1]) %in% rownames(uPExGOBP_GTEX),][1:10,],pageRank_GO_GTEX[paste(pageRank_GO_GTEX[,1]) %in% rownames(uPExGOMF_GTEX),][1:10,], pageRank_GO_GTEX[paste(pageRank_GO_GTEX[,1]) %in% rownames(uPExGOCC_GTEX),][1:10,]), Dataset="GTEX"))
pageRank_res <- cbind(GOdb[paste(pageRank_res[,1]),1:3],pageRank_res[,-1])
colnames(pageRank_res)[c(1:4)] <- c("ID","DB","Term","PRscore")
pageRank_res <- rbind(pageRank_res,data.frame(ID=pageRank_Dis_CCLE[1:10,1], DB="DISEASE", Term="", PRscore=pageRank_Dis_CCLE[1:10,2],Dataset="CCLE"),data.frame(ID=pageRank_Dis_TCGA[1:10,1], DB="DISEASE", Term="", PRscore=pageRank_Dis_TCGA[1:10,2],Dataset="TCGA"), data.frame(ID=pageRank_Dis_GTEX[1:10,1], DB="DISEASE", Term="", PRscore=pageRank_Dis_GTEX[1:10,2],Dataset="GTEX"), data.frame(ID=pageRank_MSig_CCLE[1:10,1], DB="MSigDB", Term="", PRscore=pageRank_MSig_CCLE[1:10,2],Dataset="CCLE"),data.frame(ID=pageRank_MSig_TCGA[1:10,1], DB="MSigDB", Term="", PRscore=pageRank_MSig_TCGA[1:10,2],Dataset="TCGA"), data.frame(ID=pageRank_MSig_GTEX[1:10,1], DB="MSigDB", Term="", PRscore=pageRank_MSig_GTEX[1:10,2],Dataset="GTEX")) 
write.table(pageRank_res,"DataValidation_ELIv2_Ago20.txt", row.names=F, sep="\t", quote=F)

paper_table <- data.frame(Experiment=rep(c("CCLE","GTEX","TCGA", "TotalXDB"),each=5), DB=rep(c("GO_BP","GO_MF","GO_CC","DISEASE","MSigDB"),4))
paper_table2 <- cbind(paper_table, Mean_PRScore=c(as.numeric(sapply(paste(unique(paper_table[,1]))[-4],FUN=function(x) c(mean(pageRank_res[pageRank_res[,5]==x,4][1:10]), mean(pageRank_res[pageRank_res[,5]==x,4][11:20]), mean(pageRank_res[pageRank_res[,5]==x,4][21:30]), mean(pageRank_res[pageRank_res[,5]==x,4][31:40]), mean(pageRank_res[pageRank_res[,5]==x,4][41:50])))),c(mean(c(pageRank_res[pageRank_res[,5]=="CCLE",4][1:10], pageRank_res[pageRank_res[,5]=="GTEX",4][1:10], pageRank_res[pageRank_res[,5]=="TCGA",4][1:10])), mean(c(pageRank_res[pageRank_res[,5]=="CCLE",4][11:20], pageRank_res[pageRank_res[,5]=="GTEX",4][11:20], pageRank_res[pageRank_res[,5]=="TCGA",4][11:20])), mean(c(pageRank_res[pageRank_res[,5]=="CCLE",4][21:30], pageRank_res[pageRank_res[,5]=="GTEX",4][21:30], pageRank_res[pageRank_res[,5]=="TCGA",4][21:30])), mean(c(pageRank_res[pageRank_res[,5]=="CCLE",4][31:40], pageRank_res[pageRank_res[,5]=="GTEX",4][31:40], pageRank_res[pageRank_res[,5]=="TCGA",4][31:40])), mean(c(pageRank_res[pageRank_res[,5]=="CCLE",4][41:50], pageRank_res[pageRank_res[,5]=="GTEX",4][41:50], pageRank_res[pageRank_res[,5]=="TCGA",4][41:50])))), log2Mean_PRScore=c(as.numeric(sapply(paste(unique(paper_table[,1]))[-4],FUN=function(x) c(mean(log2(pageRank_res[pageRank_res[,5]==x,4][1:10])), mean(log2(pageRank_res[pageRank_res[,5]==x,4][11:20])), mean(log2(pageRank_res[pageRank_res[,5]==x,4][21:30])), mean(log2(pageRank_res[pageRank_res[,5]==x,4][31:40])),mean(log2(pageRank_res[pageRank_res[,5]==x,4][41:50]))))),c(mean(c(log2(pageRank_res[pageRank_res[,5]=="CCLE",4][1:10]), log2(pageRank_res[pageRank_res[,5]=="GTEX",4][1:10]), log2(pageRank_res[pageRank_res[,5]=="TCGA",4][1:10]))), mean(c(log2(pageRank_res[pageRank_res[,5]=="CCLE",4][11:20]), log2(pageRank_res[pageRank_res[,5]=="GTEX",4][11:20]), log2(pageRank_res[pageRank_res[,5]=="TCGA",4][11:20]))), mean(c(log2(pageRank_res[pageRank_res[,5]=="CCLE",4][21:30]), log2(pageRank_res[pageRank_res[,5]=="GTEX",4][21:30]), log2(pageRank_res[pageRank_res[,5]=="TCGA",4][21:30]))), mean(c(log2(pageRank_res[pageRank_res[,5]=="CCLE",4][31:40]), log2(pageRank_res[pageRank_res[,5]=="GTEX",4][41:50]), log2(pageRank_res[pageRank_res[,5]=="TCGA",4][41:50]))), mean(c(log2(pageRank_res[pageRank_res[,5]=="CCLE",4][41:50]), log2(pageRank_res[pageRank_res[,5]=="GTEX",4][41:50]), log2(pageRank_res[pageRank_res[,5]=="TCGA",4][41:50]))))))
paper_table2 <- cbind(paper_table2, uPE1_withInferences=c(sum(colSums(uPExGOBP_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10],])>0),sum(colSums(uPExGOMF_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20],])>0), sum(colSums(uPExGOCC_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30],])>0),sum(colSums(uPExDisease_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40],])>0), sum(colSums(uPExMsigDB_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50],])>0),
	sum(colSums(uPExGOBP_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10],])>0),sum(colSums(uPExGOMF_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20],])>0), sum(colSums(uPExGOCC_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30],])>0),sum(colSums(uPExDisease_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40],])>0),sum(colSums(uPExMsigDB_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50],])>0),
	sum(colSums(uPExGOBP_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10],])>0),sum(colSums(uPExGOMF_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20],])>0), sum(colSums(uPExGOCC_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30],])>0),sum(colSums(uPExDisease_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40],])>0),sum(colSums(uPExMsigDB_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50],])>0)
, length(unique(c(names(colSums(uPExGOBP_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10],]))[colSums(uPExGOBP_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10],])>0], names(colSums(uPExGOBP_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10],]))[colSums(uPExGOBP_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10],])>0],names(colSums(uPExGOBP_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10],]))[colSums(uPExGOBP_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10],])>0])))
, length(unique(c(names(colSums(uPExGOMF_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20],]))[colSums(uPExGOMF_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20],])>0], names(colSums(uPExGOMF_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20],]))[colSums(uPExGOMF_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20],])>0],names(colSums(uPExGOMF_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20],]))[colSums(uPExGOMF_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20],])>0])))
,length(unique(c(names(colSums(uPExGOCC_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30],]))[colSums(uPExGOCC_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30],])>0], names(colSums(uPExGOCC_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30],]))[colSums(uPExGOCC_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30],])>0],names(colSums(uPExGOCC_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30],]))[colSums(uPExGOCC_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30],])>0]))),
length(unique(c(names(colSums(uPExDisease_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40],]))[colSums(uPExDisease_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40],])>0], names(colSums(uPExDisease_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40],]))[colSums(uPExDisease_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40],])>0],names(colSums(uPExDisease_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40],]))[colSums(uPExDisease_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40],])>0]))),
length(unique(c(names(colSums(uPExMsigDB_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50],]))[colSums(uPExMsigDB_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50],])>0], names(colSums(uPExMsigDB_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50],]))[colSums(uPExMsigDB_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50],])>0],names(colSums(uPExMsigDB_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50],]))[colSums(uPExMsigDB_TCGA[paste(pageRank_res[pageRank_res[,4]=="TCGA",1])[41:50],])>0])))))

paper_table2 <- cbind(paper_table2, Annotated_cPE1=c(sum(rowSums(CCLE_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10])])>0),sum(rowSums(CCLE_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20])])>0),sum(rowSums(CCLE_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30])])>0),sum(rowSums(CCLE_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40])])>0),sum(rowSums(CCLE_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50]])>0),

	sum(rowSums(GTEX_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10])])>0),sum(rowSums(GTEX_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20])])>0),sum(rowSums(GTEX_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30])])>0),sum(rowSums(GTEX_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40])])>0),sum(rowSums(GTEX_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50]])>0),

	sum(rowSums(TCGA_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10])])>0),sum(rowSums(TCGA_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20])])>0),sum(rowSums(TCGA_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30])])>0),sum(rowSums(TCGA_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40])])>0),sum(rowSums(TCGA_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50]])>0),

length(unique(c(names(rowSums(CCLE_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10])]))[rowSums(CCLE_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10])])>0], names(rowSums(GTEX_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10])]))[rowSums(GTEX_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10])])>0],names(rowSums(TCGA_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10])]))[rowSums(TCGA_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10])])>0]))),
length(unique(c(names(rowSums(CCLE_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20])]))[rowSums(CCLE_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20])])>0], names(rowSums(GTEX_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20])]))[rowSums(GTEX_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20])])>0],names(rowSums(TCGA_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20])]))[rowSums(TCGA_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20])])>0]))),
length(unique(c(names(rowSums(CCLE_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30])]))[rowSums(CCLE_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30])])>0], names(rowSums(GTEX_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30])]))[rowSums(GTEX_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30])])>0],names(rowSums(TCGA_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30])]))[rowSums(TCGA_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30])])>0]))),
length(unique(c(names(rowSums(CCLE_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40])]))[rowSums(CCLE_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40])])>0], names(rowSums(GTEX_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40])]))[rowSums(GTEX_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40])])>0],names(rowSums(TCGA_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40])]))[rowSums(TCGA_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40])])>0]))),
length(unique(c(names(rowSums(CCLE_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50]]))[rowSums(CCLE_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50]])>0], names(rowSums(GTEX_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50]]))[rowSums(GTEX_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50]])>0],names(rowSums(TCGA_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50]]))[rowSums(TCGA_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50]])>0])))))
paper_table2[,1] <- factor(paper_table2[,1], levels=c("CCLE","GTEX","TCGA","TotalXDB","Total"))
paper_table2[,2] <- factor(paper_table2[,2], levels=c("GO_BP","GO_MF","GO_CC","MSigDB","DISEASE","All_DB"))

paper_table2 <- rbind(paper_table2, c("Total","All_DB", mean(pageRank_res[,4]), mean(log2(pageRank_res[,4])), 
length(unique(c(names(colSums(uPExGOBP_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10],]))[colSums(uPExGOBP_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10],])>0], names(colSums(uPExGOBP_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10],]))[colSums(uPExGOBP_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10],])>0],names(colSums(uPExGOBP_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10],]))[colSums(uPExGOBP_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10],])>0], names(colSums(uPExGOMF_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20],]))[colSums(uPExGOMF_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20],])>0], names(colSums(uPExGOMF_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20],]))[colSums(uPExGOMF_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20],])>0],names(colSums(uPExGOMF_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20],]))[colSums(uPExGOMF_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20],])>0],names(colSums(uPExGOCC_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30],]))[colSums(uPExGOCC_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30],])>0], names(colSums(uPExGOCC_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30],]))[colSums(uPExGOCC_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30],])>0],names(colSums(uPExGOCC_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30],]))[colSums(uPExGOCC_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30],])>0],names(colSums(uPExDisease_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40],]))[colSums(uPExDisease_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40],])>0], names(colSums(uPExDisease_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40],]))[colSums(uPExDisease_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40],])>0],names(colSums(uPExDisease_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40],]))[colSums(uPExDisease_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40],])>0],names(colSums(uPExMsigDB_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50],]))[colSums(uPExMsigDB_CCLE[paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50],])>0], names(colSums(uPExMsigDB_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50],]))[colSums(uPExMsigDB_GTEX[paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50],])>0],names(colSums(uPExMsigDB_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50],]))[colSums(uPExMsigDB_TCGA[paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50],])>0]))),
length(unique(c(names(rowSums(CCLE_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10])]))[rowSums(CCLE_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[1:10])])>0], names(rowSums(GTEX_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10])]))[rowSums(GTEX_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[1:10])])>0],names(rowSums(TCGA_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10])]))[rowSums(TCGA_geneXGOBP[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[1:10])])>0],names(rowSums(CCLE_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20])]))[rowSums(CCLE_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[11:20])])>0], names(rowSums(GTEX_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20])]))[rowSums(GTEX_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[11:20])])>0],names(rowSums(TCGA_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20])]))[rowSums(TCGA_geneXGOMF[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[11:20])])>0],names(rowSums(CCLE_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30])]))[rowSums(CCLE_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[21:30])])>0], names(rowSums(GTEX_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30])]))[rowSums(GTEX_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[21:30])])>0],names(rowSums(TCGA_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30])]))[rowSums(TCGA_geneXGOCC[,gsub(":",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[21:30])])>0],names(rowSums(CCLE_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40])]))[rowSums(CCLE_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[31:40])])>0], names(rowSums(GTEX_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40])]))[rowSums(GTEX_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[31:40])])>0],names(rowSums(TCGA_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40])]))[rowSums(TCGA_geneXDisease[,gsub(" ",".",paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[31:40])])>0],names(rowSums(CCLE_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50]]))[rowSums(CCLE_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="CCLE",1])[41:50]])>0], names(rowSums(GTEX_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50]]))[rowSums(GTEX_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="GTEX",1])[41:50]])>0],names(rowSums(TCGA_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50]]))[rowSums(TCGA_geneXMSigDB[,paste(pageRank_res[pageRank_res[,5]=="TCGA",1])[41:50]])>0]))) ))
write.table(paper_table2, file="UPEFinder_SummaryTable.txt", row.names=F, sep="\t", quote=F)

pdf("PE1Expression_VennDiagrams_Ago20.pdf",14,7)
par(mfrow=c(1,2))
compare3List(colnames(uPExGOBP_CCLE), colnames(uPExGOBP_GTEX), colnames(uPExGOBP_TCGA), "CCLE","GTEx","TCGA", "Expressed uPE1 (Total unknown proteins=1813)")
compare3List(rownames(CCLE_geneXGOCC), rownames(GTEX_geneXGOCC), rownames(TCGA_geneXGOCC), "CCLE","GTEx","TCGA", "Expressed cPE1 (Total cPE1=16620)")
dev.off()

pdf("uPE1withFunctionsByExp_VennDiagrams_Ago20.pdf",12,4)
par(mfrow=c(1,3))
compare5List(colnames(uPExGOBP_CCLE)[colSums(uPExGOBP_CCLE)>0], colnames(uPExGOMF_CCLE)[colSums(uPExGOMF_CCLE)>0], colnames(uPExGOCC_CCLE)[colSums(uPExGOCC_CCLE)>0], colnames(uPExMsigDB_CCLE)[colSums(uPExMsigDB_CCLE)>0], colnames(uPExDisease_CCLE)[colSums(uPExDisease_CCLE)>0],"GO_BP","GO_MF","GO_CC","MSigdB", "DISEASE", "uPE1 with infered functions in CCLE")
compare5List(colnames(uPExGOBP_GTEX)[colSums(uPExGOBP_GTEX)>0], colnames(uPExGOMF_GTEX)[colSums(uPExGOMF_GTEX)>0], colnames(uPExGOCC_GTEX)[colSums(uPExGOCC_GTEX)>0], colnames(uPExMsigDB_GTEX)[colSums(uPExMsigDB_GTEX)>0], colnames(uPExDisease_GTEX)[colSums(uPExDisease_GTEX)>0],"GO_BP","GO_MF", "GO_CC","MSigdB", "DISEASE", "uPE1 with infered functions in GTEx")
compare5List(colnames(uPExGOBP_TCGA)[colSums(uPExGOBP_TCGA)>0], colnames(uPExGOMF_TCGA)[colSums(uPExGOMF_TCGA)>0], colnames(uPExGOCC_TCGA)[colSums(uPExGOCC_TCGA)>0], colnames(uPExMsigDB_TCGA)[colSums(uPExMsigDB_TCGA)>0], colnames(uPExDisease_TCGA)[colSums(uPExDisease_TCGA)>0],"GO_BP","GO_MF","GO_CC","MSigdB", "DISEASE", "uPE1 with infered functions in TCGA")
dev.off()

pdf("uPE1withFunctionsByDB_VennDiagrams_Ago20.pdf",12,12)
par(mfrow=c(2,2))
compare3List(colnames(uPExGOBP_CCLE)[colSums(uPExGOBP_CCLE)>0], colnames(uPExGOBP_GTEX)[colSums(uPExGOBP_GTEX)>0], colnames(uPExGOBP_TCGA)[colSums(uPExGOBP_TCGA)>0], "CCLE","GTEx","TCGA", "uPE1 with infered functions with GO_BP")
compare3List(colnames(uPExGOMF_CCLE)[colSums(uPExGOMF_CCLE)>0], colnames(uPExGOMF_GTEX)[colSums(uPExGOMF_GTEX)>0], colnames(uPExGOMF_TCGA)[colSums(uPExGOMF_TCGA)>0], "CCLE","GTEx","TCGA", "uPE1 with infered functions with GO_MF")
compare3List(colnames(uPExGOCC_CCLE)[colSums(uPExGOCC_CCLE)>0], colnames(uPExGOCC_GTEX)[colSums(uPExGOCC_GTEX)>0], colnames(uPExGOCC_TCGA)[colSums(uPExGOCC_TCGA)>0], "CCLE","GTEx","TCGA", "uPE1 with infered functions with GO_CC")
compare3List(colnames(uPExMsigDB_CCLE)[colSums(uPExMsigDB_CCLE)>0], colnames(uPExMsigDB_GTEX)[colSums(uPExMsigDB_GTEX)>0], colnames(uPExMsigDB_TCGA)[colSums(uPExMsigDB_TCGA)>0], "CCLE","GTEx","TCGA", "uPE1 with infered functions with MSigDB")
compare3List(colnames(uPExDisease_CCLE)[colSums(uPExDisease_CCLE)>0], colnames(uPExDisease_GTEX)[colSums(uPExDisease_GTEX)>0], colnames(uPExDisease_TCGA)[colSums(uPExDisease_TCGA)>0], "CCLE","GTEx","TCGA", "uPE1 with infered functions with DISEASE")
dev.off()

pdf("uPE1withFunctionsTotal_VennDiagrams_Ago20.pdf",8,8)
compare5List(unique(c(colnames(uPExGOBP_CCLE)[colSums(uPExGOBP_CCLE)>0], colnames(uPExGOBP_GTEX)[colSums(uPExGOBP_GTEX)>0], colnames(uPExGOBP_TCGA)[colSums(uPExGOBP_TCGA)>0])), unique(c(colnames(uPExGOMF_CCLE)[colSums(uPExGOMF_CCLE)>0], colnames(uPExGOMF_GTEX)[colSums(uPExGOMF_GTEX)>0], colnames(uPExGOMF_TCGA)[colSums(uPExGOMF_TCGA)>0])), unique(c(colnames(uPExGOCC_CCLE)[colSums(uPExGOCC_CCLE)>0], colnames(uPExGOCC_GTEX)[colSums(uPExGOCC_GTEX)>0], colnames(uPExGOCC_TCGA)[colSums(uPExGOCC_TCGA)>0])), unique(c(colnames(uPExMsigDB_CCLE)[colSums(uPExMsigDB_CCLE)>0], colnames(uPExMsigDB_GTEX)[colSums(uPExMsigDB_GTEX)>0], colnames(uPExMsigDB_TCGA)[colSums(uPExMsigDB_TCGA)>0])), unique(c(colnames(uPExDisease_CCLE)[colSums(uPExDisease_CCLE)>0], colnames(uPExDisease_GTEX)[colSums(uPExDisease_GTEX)>0], colnames(uPExDisease_TCGA)[colSums(uPExDisease_TCGA)>0])),"GO_BP", "GO_MF","GO_CC","MSigdB", "DISEASE", "Total uPE1 with infered functions by DB")

compare3List(unique(c(colnames(uPExGOBP_CCLE)[colSums(uPExGOBP_CCLE)>0], colnames(uPExGOMF_CCLE)[colSums(uPExGOMF_CCLE)>0],colnames(uPExGOCC_CCLE)[colSums(uPExGOCC_CCLE)>0],colnames(uPExMsigDB_CCLE)[colSums(uPExMsigDB_CCLE)>0],colnames(uPExDisease_CCLE)[colSums(uPExDisease_CCLE)>0])), unique(c(colnames(uPExGOBP_GTEX)[colSums(uPExGOBP_GTEX)>0], colnames(uPExGOMF_GTEX)[colSums(uPExGOMF_GTEX)>0], colnames(uPExGOCC_GTEX)[colSums(uPExGOCC_GTEX)>0],colnames(uPExMsigDB_GTEX)[colSums(uPExMsigDB_GTEX)>0],colnames(uPExDisease_GTEX)[colSums(uPExDisease_GTEX)>0])), unique(c(colnames(uPExGOBP_TCGA)[colSums(uPExGOBP_TCGA)>0],colnames(uPExGOCC_TCGA)[colSums(uPExGOCC_TCGA)>0],colnames(uPExMsigDB_TCGA)[colSums(uPExMsigDB_TCGA)>0],colnames(uPExDisease_TCGA)[colSums(uPExDisease_TCGA)>0])), "CCLE","GTEx","TCGA", "uPE1 with infered functions by dataset")

dev.off()


## uPE DATA
uPExGOBP_TCGA <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_enriched_GOBP_filtered.txt", header=T, sep="\t")
rownames(uPExGOBP_TCGA) <- uPExGOBP_TCGA[,1]
uPExGOBP_TCGA <- uPExGOBP_TCGA[,-1]
uPExGOCC_TCGA <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_enriched_GOCC_filtered.txt", header=T, sep="\t")
rownames(uPExGOCC_TCGA) <- uPExGOCC_TCGA[,1]
uPExGOCC_TCGA <- uPExGOCC_TCGA[,-1]
uPExDisease_TCGA <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_enriched_Disease_filtered.txt", header=T, sep="\t")
rownames(uPExDisease_TCGA) <- uPExDisease_TCGA[,1]
uPExDisease_TCGA <- uPExDisease_TCGA[,-1]
uPExMsigDB_TCGA <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_enriched_MSigDB.txt", header=T, sep="\t")
rownames(uPExMsigDB_TCGA) <- uPExMsigDB_TCGA[,1]
uPExMsigDB_TCGA <- uPExMsigDB_TCGA[,-1]
uPExGOMF_TCGA <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_enriched_GOMF_filtered.txt", header=T, sep="\t")
rownames(uPExGOMF_TCGA) <- uPExGOMF_TCGA[,1]
uPExGOMF_TCGA <- uPExGOMF_TCGA[,-1]

uPExGOBP_CCLE <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_enriched_GOBP_filtered.txt", header=T, sep="\t")
rownames(uPExGOBP_CCLE) <- uPExGOBP_CCLE[,1]
uPExGOBP_CCLE <- uPExGOBP_CCLE[,-1]
uPExGOCC_CCLE <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_enriched_GOCC_filtered.txt", header=T, sep="\t")
rownames(uPExGOCC_CCLE) <- uPExGOCC_CCLE[,1]
uPExGOCC_CCLE <- uPExGOCC_CCLE[,-1]
uPExDisease_CCLE <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_enriched_Disease_filtered.txt", header=T, sep="\t")
rownames(uPExDisease_CCLE) <- uPExDisease_CCLE[,1]
uPExDisease_CCLE <- uPExDisease_CCLE[,-1]
uPExMsigDB_CCLE <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_enriched_MSigDB.txt", header=T, sep="\t")
rownames(uPExMsigDB_CCLE) <- uPExMsigDB_CCLE[,1]
uPExMsigDB_CCLE <- uPExMsigDB_CCLE[,-1]
uPExGOMF_CCLE <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_enriched_GOMF_filtered.txt", header=T, sep="\t")
rownames(uPExGOMF_CCLE) <- uPExGOMF_CCLE[,1]
uPExGOMF_CCLE <- uPExGOMF_CCLE[,-1]

uPExGOBP_GTEX <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_enriched_GOBP_filtered.txt", header=T, sep="\t")
rownames(uPExGOBP_GTEX) <- uPExGOBP_GTEX[,1]
uPExGOBP_GTEX <- uPExGOBP_GTEX[,-1]
uPExGOCC_GTEX <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_enriched_GOCC_filtered.txt", header=T, sep="\t")
rownames(uPExGOCC_GTEX) <- uPExGOCC_GTEX[,1]
uPExGOCC_GTEX <- uPExGOCC_GTEX[,-1]
uPExDisease_GTEX <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_enriched_Disease_filtered.txt", header=T, sep="\t")
rownames(uPExDisease_GTEX) <- uPExDisease_GTEX[,1]
uPExDisease_GTEX <- uPExDisease_GTEX[,-1]
uPExMsigDB_GTEX <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_enriched_MSigDB.txt", header=T, sep="\t")
rownames(uPExMsigDB_GTEX) <- uPExMsigDB_GTEX[,1]
uPExMsigDB_GTEX <- uPExMsigDB_GTEX[,-1]
uPExGOMF_GTEX <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_enriched_GOMF_filtered.txt", header=T, sep="\t")
rownames(uPExGOMF_GTEX) <- uPExGOMF_GTEX[,1]
uPExGOMF_GTEX <- uPExGOMF_GTEX[,-1]

### Paper figures
ProteinLinks_ggplot2 <- rbind(data.frame(Protein=unlist(sapply(colnames(uPExGOBP_TCGA), FUN=function(x) rep(x, sum(uPExGOBP_TCGA[,x])))),ProteinClass="uPE",Experiment="TCGA",DB="GO_BP",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExGOMF_TCGA), FUN=function(x) rep(x, sum(uPExGOMF_TCGA[,x])))),ProteinClass="uPE",Experiment="TCGA",DB="GO_MF",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExGOCC_TCGA), FUN=function(x) rep(x, sum(uPExGOCC_TCGA[,x])))),ProteinClass="uPE",Experiment="TCGA",DB="GO_CC",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExDisease_TCGA), FUN=function(x) rep(x, sum(uPExDisease_TCGA[,x])))),ProteinClass="uPE",Experiment="TCGA",DB="Disease",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExMsigDB_TCGA), FUN=function(x) rep(x, sum(uPExMsigDB_TCGA[,x])))),ProteinClass="uPE",Experiment="TCGA",DB="MSigDB",Link="Out"))
ProteinLinks_ggplot2 <- rbind(ProteinLinks_ggplot2,data.frame(Protein=unlist(sapply(colnames(uPExGOBP_CCLE), FUN=function(x) rep(x, sum(uPExGOBP_CCLE[,x])))),ProteinClass="uPE",Experiment="CCLE",DB="GO_BP",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExGOMF_CCLE), FUN=function(x) rep(x, sum(uPExGOMF_CCLE[,x])))),ProteinClass="uPE",Experiment="CCLE",DB="GO_MF",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExGOCC_CCLE), FUN=function(x) rep(x, sum(uPExGOCC_CCLE[,x])))),ProteinClass="uPE",Experiment="CCLE",DB="GO_CC",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExDisease_CCLE), FUN=function(x) rep(x, sum(uPExDisease_CCLE[,x])))),ProteinClass="uPE",Experiment="CCLE",DB="Disease",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExMsigDB_CCLE), FUN=function(x) rep(x, sum(uPExMsigDB_CCLE[,x])))),ProteinClass="uPE",Experiment="CCLE",DB="MSigDB",Link="Out"))
ProteinLinks_ggplot2 <- rbind(ProteinLinks_ggplot2,data.frame(Protein=unlist(sapply(colnames(uPExGOBP_GTEX), FUN=function(x) rep(x, sum(uPExGOBP_GTEX[,x])))),ProteinClass="uPE",Experiment="GTEX",DB="GO_BP",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExGOMF_GTEX), FUN=function(x) rep(x, sum(uPExGOMF_GTEX[,x])))),ProteinClass="uPE",Experiment="GTEX",DB="GO_MF",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExGOCC_GTEX), FUN=function(x) rep(x, sum(uPExGOCC_GTEX[,x])))),ProteinClass="uPE",Experiment="GTEX",DB="GO_CC",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExDisease_GTEX), FUN=function(x) rep(x, sum(uPExDisease_GTEX[,x])))),ProteinClass="uPE",Experiment="GTEX",DB="Disease",Link="Out"), data.frame(Protein=unlist(sapply(colnames(uPExMsigDB_GTEX), FUN=function(x) rep(x, sum(uPExMsigDB_GTEX[,x])))),ProteinClass="uPE",Experiment="GTEX",DB="MSigDB",Link="Out"))
ProteinLinks_ggplot2 <- rbind(ProteinLinks_ggplot2,data.frame(Protein=unlist(sapply(colnames(TCGA_uPExPE1), FUN=function(x) rep(x, sum(TCGA_uPExPE1[,x])))),ProteinClass="uPE",Experiment="TCGA",DB="Cor",Link="Out"), data.frame(Protein=unlist(sapply(colnames(CCLE_uPExPE1), FUN=function(x) rep(x, sum(CCLE_uPExPE1[,x])))),ProteinClass="uPE",Experiment="CCLE",DB="Cor",Link="Out"), data.frame(Protein=unlist(sapply(colnames(GTEX_uPExPE1), FUN=function(x) rep(x, sum(GTEX_uPExPE1[,x])))),ProteinClass="uPE",Experiment="GTEX",DB="Cor",Link="Out"))


ProteinLinks2_ggplot2 <- rbind(data.frame(Protein=unlist(sapply(rownames(TCGA_geneXGOBP), FUN=function(x) rep(x, sum(TCGA_geneXGOBP[x,])))),ProteinClass="cPE1",Experiment="TCGA",DB="GO_BP",Link="Out"), data.frame(Protein=unlist(sapply(rownames(TCGA_geneXGOMF), FUN=function(x) rep(x, sum(TCGA_geneXGOMF[x,])))),ProteinClass="cPE1",Experiment="TCGA",DB="GO_MF",Link="Out"), data.frame(Protein=unlist(sapply(rownames(TCGA_geneXGOCC), FUN=function(x) rep(x, sum(TCGA_geneXGOCC[x,])))),ProteinClass="cPE1",Experiment="TCGA",DB="GO_CC",Link="Out"), data.frame(Protein=unlist(sapply(rownames(TCGA_geneXDisease), FUN=function(x) rep(x, sum(TCGA_geneXDisease[x,])))),ProteinClass="cPE1",Experiment="TCGA",DB="Disease",Link="Out"), data.frame(Protein=unlist(sapply(rownames(TCGA_geneXMSigDB), FUN=function(x) rep(x, sum(TCGA_geneXMSigDB[x,])))),ProteinClass="cPE1",Experiment="TCGA",DB="MSigDB",Link="Out"))
ProteinLinks2_ggplot2 <- rbind(ProteinLinks2_ggplot2,data.frame(Protein=unlist(sapply(rownames(CCLE_geneXGOBP), FUN=function(x) rep(x, sum(CCLE_geneXGOBP[x,])))),ProteinClass="cPE1",Experiment="CCLE",DB="GO_BP",Link="Out"), data.frame(Protein=unlist(sapply(rownames(CCLE_geneXGOMF), FUN=function(x) rep(x, sum(CCLE_geneXGOMF[x,])))),ProteinClass="cPE1",Experiment="CCLE",DB="GO_MF",Link="Out"), data.frame(Protein=unlist(sapply(rownames(CCLE_geneXGOCC), FUN=function(x) rep(x, sum(CCLE_geneXGOCC[x,])))),ProteinClass="cPE1",Experiment="CCLE",DB="GO_CC",Link="Out"), data.frame(Protein=unlist(sapply(rownames(CCLE_geneXDisease), FUN=function(x) rep(x, sum(CCLE_geneXDisease[x,])))),ProteinClass="cPE1",Experiment="CCLE",DB="Disease",Link="Out"), data.frame(Protein=unlist(sapply(rownames(CCLE_geneXMSigDB), FUN=function(x) rep(x, sum(CCLE_geneXMSigDB[x,])))),ProteinClass="cPE1",Experiment="CCLE",DB="MSigDB",Link="Out"))
ProteinLinks2_ggplot2 <- rbind(ProteinLinks2_ggplot2,data.frame(Protein=unlist(sapply(rownames(GTEX_geneXGOBP), FUN=function(x) rep(x, sum(GTEX_geneXGOBP[x,])))),ProteinClass="cPE1",Experiment="GTEX",DB="GO_BP",Link="Out"), data.frame(Protein=unlist(sapply(rownames(GTEX_geneXGOMF), FUN=function(x) rep(x, sum(GTEX_geneXGOMF[x,])))),ProteinClass="cPE1",Experiment="GTEX",DB="GO_MF",Link="Out"), data.frame(Protein=unlist(sapply(rownames(GTEX_geneXGOCC), FUN=function(x) rep(x, sum(GTEX_geneXGOCC[x,])))),ProteinClass="cPE1",Experiment="GTEX",DB="GO_CC",Link="Out"), data.frame(Protein=unlist(sapply(rownames(GTEX_geneXDisease), FUN=function(x) rep(x, sum(GTEX_geneXDisease[x,])))),ProteinClass="cPE1",Experiment="GTEX",DB="Disease",Link="Out"), data.frame(Protein=unlist(sapply(rownames(GTEX_geneXMSigDB), FUN=function(x) rep(x, sum(GTEX_geneXMSigDB[x,])))),ProteinClass="cPE1",Experiment="GTEX",DB="MSigDB",Link="Out"))
ProteinLinks2_ggplot2 <- rbind(ProteinLinks2_ggplot2,data.frame(Protein=unlist(sapply(rownames(TCGA_uPExPE1), FUN=function(x) rep(x, sum(TCGA_uPExPE1[x,])))),ProteinClass="cPE1",Experiment="TCGA",DB="Cor",Link="In"), data.frame(Protein=unlist(sapply(rownames(CCLE_uPExPE1), FUN=function(x) rep(x, sum(CCLE_uPExPE1[x,])))),ProteinClass="cPE1",Experiment="CCLE",DB="Cor",Link="In"), data.frame(Protein=unlist(sapply(rownames(GTEX_uPExPE1), FUN=function(x) rep(x, sum(GTEX_uPExPE1[x,])))),ProteinClass="cPE1",Experiment="GTEX",DB="Cor",Link="In"))
ProteinLinks2_ggplot2$Link <- factor(paste(ProteinLinks2_ggplot2$Link), levels=c("Out","In"))

DBLinks_ggplot2 <- rbind(data.frame(Protein=unlist(sapply(rownames(uPExGOBP_TCGA), FUN=function(x) rep(x, sum(uPExGOBP_TCGA[x,])))),ProteinClass="uPE",Experiment="TCGA",DB="GO_BP",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExGOMF_TCGA), FUN=function(x) rep(x, sum(uPExGOMF_TCGA[x,])))),ProteinClass="uPE",Experiment="TCGA",DB="GO_MF",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExGOCC_TCGA), FUN=function(x) rep(x, sum(uPExGOCC_TCGA[x,])))),ProteinClass="uPE",Experiment="TCGA",DB="GO_CC",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExDisease_TCGA), FUN=function(x) rep(x, sum(uPExDisease_TCGA[x,])))),ProteinClass="uPE",Experiment="TCGA",DB="Disease",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExMsigDB_TCGA), FUN=function(x) rep(x, sum(uPExMsigDB_TCGA[x,])))),ProteinClass="uPE",Experiment="TCGA",DB="MSigDB",Link="In"))
DBLinks_ggplot2 <- rbind(DBLinks_ggplot2,data.frame(Protein=unlist(sapply(rownames(uPExGOBP_CCLE), FUN=function(x) rep(x, sum(uPExGOBP_CCLE[x,])))),ProteinClass="uPE",Experiment="CCLE",DB="GO_BP",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExGOMF_CCLE), FUN=function(x) rep(x, sum(uPExGOMF_CCLE[x,])))),ProteinClass="uPE",Experiment="CCLE",DB="GO_MF",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExGOCC_CCLE), FUN=function(x) rep(x, sum(uPExGOCC_CCLE[x,])))),ProteinClass="uPE",Experiment="CCLE",DB="GO_CC",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExDisease_CCLE), FUN=function(x) rep(x, sum(uPExDisease_CCLE[x,])))),ProteinClass="uPE",Experiment="CCLE",DB="Disease",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExMsigDB_CCLE), FUN=function(x) rep(x, sum(uPExMsigDB_CCLE[x,])))),ProteinClass="uPE",Experiment="CCLE",DB="MSigDB",Link="In"))
DBLinks_ggplot2 <- rbind(DBLinks_ggplot2,data.frame(Protein=unlist(sapply(rownames(uPExGOBP_GTEX), FUN=function(x) rep(x, sum(uPExGOBP_GTEX[x,])))),ProteinClass="uPE",Experiment="GTEX",DB="GO_BP",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExGOMF_GTEX), FUN=function(x) rep(x, sum(uPExGOMF_GTEX[x,])))),ProteinClass="uPE",Experiment="GTEX",DB="GO_MF",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExGOCC_GTEX), FUN=function(x) rep(x, sum(uPExGOCC_GTEX[x,])))),ProteinClass="uPE",Experiment="GTEX",DB="GO_CC",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExDisease_GTEX), FUN=function(x) rep(x, sum(uPExDisease_GTEX[x,])))),ProteinClass="uPE",Experiment="GTEX",DB="Disease",Link="In"), data.frame(Protein=unlist(sapply(rownames(uPExMsigDB_GTEX), FUN=function(x) rep(x, sum(uPExMsigDB_GTEX[x,])))),ProteinClass="uPE",Experiment="GTEX",DB="MSigDB",Link="In"))
DBLinks_ggplot2 <- rbind(DBLinks_ggplot2,data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(TCGA_geneXGOBP), FUN=function(x) rep(x, sum(TCGA_geneXGOBP[,x]))))),ProteinClass="cPE1",Experiment="TCGA",DB="GO_BP",Link="In"), data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(TCGA_geneXGOMF), FUN=function(x) rep(x, sum(TCGA_geneXGOMF[,x]))))),ProteinClass="cPE1",Experiment="TCGA",DB="GO_MF",Link="In"), data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(TCGA_geneXGOCC), FUN=function(x) rep(x, sum(TCGA_geneXGOCC[,x]))))),ProteinClass="cPE1",Experiment="TCGA",DB="GO_CC",Link="In"), data.frame(Protein=gsub("[.]"," ",unlist(sapply(colnames(TCGA_geneXDisease), FUN=function(x) rep(x, sum(TCGA_geneXDisease[,x]))))),ProteinClass="cPE1",Experiment="TCGA",DB="Disease",Link="In"), data.frame(Protein=unlist(sapply(colnames(TCGA_geneXMSigDB), FUN=function(x) rep(x, sum(TCGA_geneXMSigDB[,x])))),ProteinClass="cPE1",Experiment="TCGA",DB="MSigDB",Link="In"))
DBLinks_ggplot2 <- rbind(DBLinks_ggplot2,data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(CCLE_geneXGOBP), FUN=function(x) rep(x, sum(CCLE_geneXGOBP[,x]))))),ProteinClass="cPE1",Experiment="CCLE",DB="GO_BP",Link="In"), data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(CCLE_geneXGOMF), FUN=function(x) rep(x, sum(CCLE_geneXGOMF[,x]))))), ProteinClass="cPE1",Experiment="CCLE",DB="GO_MF",Link="In"), data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(CCLE_geneXGOCC), FUN=function(x) rep(x, sum(CCLE_geneXGOCC[,x]))))),ProteinClass="cPE1",Experiment="CCLE",DB="GO_CC",Link="In"), data.frame(Protein=gsub("[.]"," ",unlist(sapply(colnames(CCLE_geneXDisease), FUN=function(x) rep(x, sum(CCLE_geneXDisease[,x]))))),ProteinClass="cPE1",Experiment="CCLE",DB="Disease",Link="In"), data.frame(Protein=unlist(sapply(colnames(CCLE_geneXMSigDB), FUN=function(x) rep(x, sum(CCLE_geneXMSigDB[,x])))),ProteinClass="cPE1",Experiment="CCLE",DB="MSigDB",Link="In"))
DBLinks_ggplot2 <- rbind(DBLinks_ggplot2,data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(GTEX_geneXGOBP), FUN=function(x) rep(x, sum(GTEX_geneXGOBP[,x]))))),ProteinClass="cPE1",Experiment="GTEX",DB="GO_BP",Link="In"), data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(GTEX_geneXGOMF), FUN=function(x) rep(x, sum(GTEX_geneXGOMF[,x]))))),ProteinClass="cPE1",Experiment="GTEX",DB="GO_MF",Link="In"), data.frame(Protein=gsub("[.]",":",unlist(sapply(colnames(GTEX_geneXGOCC), FUN=function(x) rep(x, sum(GTEX_geneXGOCC[,x]))))),ProteinClass="cPE1",Experiment="GTEX",DB="GO_CC",Link="In"), data.frame(Protein=gsub("[.]"," ",unlist(sapply(colnames(GTEX_geneXDisease), FUN=function(x) rep(x, sum(GTEX_geneXDisease[,x]))))),ProteinClass="cPE1",Experiment="GTEX",DB="Disease",Link="In"), data.frame(Protein=unlist(sapply(colnames(GTEX_geneXMSigDB), FUN=function(x) rep(x, sum(GTEX_geneXMSigDB[,x])))),ProteinClass="cPE1",Experiment="GTEX",DB="MSigDB",Link="In"))

p <- ggplot(ProteinLinks_ggplot2, aes(x=Protein,color=Link, fill=Link)) + geom_bar(stat="count") + xlab("uPEs") + facet_grid(Experiment ~ DB, scale="free_x") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_colour_manual(values=c("#FC766AFF","#5884B1FF")) + scale_fill_manual(values=c("#FC766AFF","#5884B1FF"))
pdf("Links_uPEs_Ago20.pdf")
p
dev.off()
# density with categorical variable
p <- ggplot(ProteinLinks_ggplot2, aes(x=Protein,color=Link, fill=Link)) + geom_bar(stat="count", aes(y = ..prop.. , group=1)) + xlab("uPEs") + facet_grid(Experiment ~ DB, scale="free_x") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_colour_manual(values=c("#FC766AFF","#5884B1FF")) + scale_fill_manual(values=c("#FC766AFF","#5884B1FF")) 
pdf("Links_density_uPEs_Ago20.pdf")
p
dev.off()

p <- ggplot(ProteinLinks2_ggplot2, aes(x=Protein,color=Link, fill=Link)) + geom_bar(stat="count") + xlab("cPE1s") + facet_grid(Experiment ~ DB, scale="free_x") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_colour_manual(values=c("#FC766AFF","#5884B1FF")) + scale_fill_manual(values=c("#FC766AFF","#5884B1FF"))
pdf("Links_cPE1_Ago20.pdf")
p
dev.off()
# density with categorical variable
p <- ggplot(ProteinLinks2_ggplot2, aes(x=Protein,color=Link, fill=Link)) + geom_bar(stat="count", aes(y = ..prop.. , group=1)) + xlab("cPE1s") + facet_grid(Experiment ~ DB, scale="free_x") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_colour_manual(values=c("#FC766AFF","#5884B1FF")) + scale_fill_manual(values=c("#FC766AFF","#5884B1FF"))
pdf("Links_density_cPE1_Ago20.pdf")
p
dev.off()

p <- ggplot(DBLinks_ggplot2, aes(x=Protein,color=Link, fill=Link)) + geom_bar(stat="count") + xlab("DB_IDs") + facet_grid(Experiment ~ DB, scale="free_x") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_colour_manual(values=c("#5884B1FF","#FC766AFF")) + scale_fill_manual(values=c("#5884B1FF","#FC766AFF"))
pdf("Links_DB_Ago20.pdf")
p
dev.off()
# density with categorical variable
p <- ggplot(DBLinks_ggplot2, aes(x=Protein,color=Link, fill=Link)) + geom_bar(stat="count", aes(y = ..prop.. , group=1)) + xlab("DB_IDs") + facet_grid(Experiment ~ DB, scale="free_x") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_colour_manual(values=c("#5884B1FF","#FC766AFF")) + scale_fill_manual(values=c("#5884B1FF","#FC766AFF"))
pdf("Links_density_DB_Ago20.pdf")
p
dev.off()

# Comprobación con histogram/density y variable cuantitativa
a <- t(sapply(paste(unique(ProteinLinks_ggplot2[ProteinLinks_ggplot2$DB=="Cor" & ProteinLinks_ggplot2$Experiment=="TCGA", "Protein"])), FUN=function(x) c(Protein=x, Link=sum(ProteinLinks_ggplot2$DB=="Cor" & ProteinLinks_ggplot2$Experiment=="TCGA" & ProteinLinks_ggplot2$Protein==x))))
a_df <- as.data.frame(a)
a_df$Link <- as.numeric(paste(a[,2]))
p <- ggplot(a_df, aes(x=Protein, y=Link)) + geom_histogram(stat="identity") + geom_density(aes(alpha=0.2)) + xlab("uPEs") + ylab("Correlated cPE1s (TCGA)") + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
pdf("Links_uPEs_CorTCGA.pdf")
p
dev.off()

########## TCGA ##############
# Membrane (320) metabolic process(331)
colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0016020",]==1]
colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1]

# Diabetes mellitus (132) Gastrointestinal system disease (228) Neurodegenerative disease (136) Autoimmune hypersensitivity disease (135)
colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1]
colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1]
colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Neurodegenerative disease",]==1]
colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Autoimmune hypersensitivity disease",]==1]
# CHEN_METABOLIC_SYNDROM_NETWORK (160) BLALOCK_ALZHEIMERS_DISEASE_UP (133) MARSON_BOUND_BY_FOXP3_UNSTIMULATED (225)
colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1]
colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1]
colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1]

## Metabolism (62) metabolismXmembrane (51) metabolismDisXmetabolism (0)
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0016020",]==1])))
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1])))

# Cambio para el paper: Organelle membrane (296) protein metabolic process (289) Diabetes mellitus (172) Gastrointestinal system disease (306) CHEN_METABOLIC_SYNDROM_NETWORK (204) MetabolismDis (77) MetabolismDisXmembrane(2) MetabolismDisXMetabolismProcesses (0)
colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0031090",]==1]
colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0019538",]==1]
colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1] # NIBAN, TTC24 ..
colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1] # NIBAN, TTC24...
colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1] # NIBAN, TTC24
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1],colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1])) # NIBAN
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0031090",]==1])))
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1])))

## Neuro (16) NeuroXmembrane (9) NeuroXmetabolism (2) NeuroXmembraneXmetabolism (0)
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0016020",]==1]))
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1]))
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA[\"Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA[\"BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], intersect(colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1], colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0016020",]==1])))

## Immune (36) ImmuneXmembrane (33) ImmuneXmetabolism (0)
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA[\"Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0016020",]==1]))
intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1]))

########## CCLE ##############
# Membrane (488) metabolism (117)
colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0016020",]==1]
colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1]
# Diabetes mellitus (154) Gastrointestinal system disease (399) Neurodegenerative disease (212) Cerebrovascular disease (137) Autoimmune hypersensitivity disease (315)
colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1]
colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1]
colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Neurodegenerative disease",]==1]
colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Cerebrovascular disease",]==1]
colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Autoimmune hypersensitivity disease",]==1]
# CHEN_METABOLIC_SYNDROM_NETWORK (426) BLALOCK_ALZHEIMERS_DISEASE_UP (262) MARSON_BOUND_BY_FOXP3_UNSTIMULATED (130) LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN (91)
colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1]
colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1]
colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN",]==1]
colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1]


## Metabolism (118) metabolismXmembrane (101) metabolismDisXmetabolism (1) metabolismDisXmetabolismXmembrane (0)
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0016020",]==1])))
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1])))
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], intersect(colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1], colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0016020",]==1]))))

## Neuro (2) NeuroXmembrane (2) NeuroXmetabolism (0)
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN",]==1],intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Cerebrovascular disease",]==1],colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0016020",]==1]))))
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN",]==1],intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Cerebrovascular disease",]==1], colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1]))))
## Immune (59) ImmuneXmembrane (50) ImmuneXmetabolism (2) ImmuneXmetabolismXmembrane (0)
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0016020",]==1]))
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1]))
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], intersect(colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1], colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0016020",]==1])))

# Cambio para el paper: Organelle membrane (107) protein metabolic process (167) Diabetes mellitus (197) Gastrointestinal system disease (500) CHEN_METABOLIC_SYNDROM_NETWORK (530) MetabolismDis (150) MetabolismDisXmembrane(19) MetabolismDisXMetabolismProcesses (48) MetabolismDisXmembraneXMetabolicProcesses (14)
colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0031090",]==1] # NIBAN y TTC24
colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0019538",]==1]
colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1] # NIBAN, TTC24...
colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1] # NIBAN, TTC24...
colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1] # NIBAN, TTC24 ...
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1])) # NIBAN
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0031090",]==1]))) # NIBAN
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0019538",]==1])))
intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], intersect(colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0019538",]==1], colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0031090",]==1]))))


########## GTEX ##############
# Membrane (508) Metabolic process (157)
colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1]
colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1]
# Diabetes mellitus (101) Gastrointestinal system disease (176) Neurodegenerative disease (357) Autoimmune hypersensitivity disease (244)
colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1]
colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1]
colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Neurodegenerative disease",]==1]
colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Autoimmune hypersensitivity disease",]==1]
# CHEN_METABOLIC_SYNDROM_NETWORK (313) BLALOCK_ALZHEIMERS_DISEASE_UP (205) MARSON_BOUND_BY_FOXP3_UNSTIMULATED (174)
colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1]
colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1]
colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1]


## Metabolism (62) metabolismXmembrane (52) metabolismDisXmetabolism (2) metabolismDisXmetabolismXmembrane (0)
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1])))
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1])))
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], intersect(colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1],colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1]))))
## Neuro (102) NeuroXmembrane (86) NeuroXmetabolism (2) NeuroXmetabolismXmembrane (0)
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1]))
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1]))
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], intersect(colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1], colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1])))
## ImmuneXmembrane (44) ImmuneXmetabolism (5) ImmuneXmetabolismXmembrane (0)
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1],colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1]))
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1],colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1]))
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1],intersect(colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1],colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1])))

# Cambio para el paper: Organelle membrane (213) protein metabolic process (209) Diabetes mellitus (144) Gastrointestinal system disease (216) CHEN_METABOLIC_SYNDROM_NETWORK (309) MetabolismDis (81) MetabolismDisXmembrane(0) MetabolismDisXMetabolismProcesses (7) MetabolismDisXmembraneXMetabolicProcesses (0)
colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0031090",]==1] 
colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0019538",]==1]
colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1] # NIBAN, TTC24.. (but not SOWAHD)
colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1] # NIBAN, TTC24..
colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1] # NIBAN, TTC24..
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1])) # NIBAN
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0031090",]==1]))) 
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0019538",]==1])))
intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], intersect(colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0019538",]==1], colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0031090",]==1]))))



#-----------
# 5 metabXmembrane
metabolism_uPE <- intersect(intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0016020",]==1]))), intersect(intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0016020",]==1]))), intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1])))))
# 0 metabXmetabolism
metabolism_uPE_2 <- intersect(intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1]))), intersect(intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1]))), intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1])))))
# 0 neuroXmembrane
Neuro_uPE <- intersect(intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0016020",]==1])),intersect(intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN",]==1],intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Cerebrovascular disease",]==1],colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0016020",]==1])))),intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1]))))
# 0 neuroXmetabolism
Neuro_uPE_2 <- intersect(intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1])),intersect(intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN",]==1],intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Cerebrovascular disease",]==1], colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1])))),intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Neurodegenerative disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["BLALOCK_ALZHEIMERS_DISEASE_UP",]==1], colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1]))))
# 2 inmmuneXmembrane
Immune_uPE <- intersect(intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], colnames(uPExGOCC_TCGA)[uPExGOCC_TCGA["GO:0016020",]==1])),intersect(intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0016020",]==1])),
	intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1],colnames(uPExGOCC_GTEX)[uPExGOCC_GTEX["GO:0016020",]==1]))))
# 0 inmmuneXmetabolism
Immune_uPE_2 <- intersect(intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], colnames(uPExGOBP_TCGA)[uPExGOBP_TCGA["GO:0008152",]==1])),intersect(intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1], colnames(uPExGOBP_CCLE)[uPExGOBP_CCLE["GO:0008152",]==1]))
	,intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Autoimmune hypersensitivity disease",]==1], intersect(colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["MARSON_BOUND_BY_FOXP3_UNSTIMULATED",]==1],colnames(uPExGOBP_GTEX)[uPExGOBP_GTEX["GO:0008152",]==1]))))

# Para el paper
# 11 metabolismDis (4 de las 5 que tenía antes (excepto SOWAHD))
metabolism_uPE <- intersect(intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1], colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1])), intersect(intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1])), intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1]))))

 intersect(intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_TCGA)[uPExDisease_TCGA["Gastrointestinal system disease",]==1], colnames(uPExMsigDB_TCGA)[uPExMsigDB_TCGA["CHEN_METABOLIC_SYNDROM_NETWORK",]==1])), intersect(intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_CCLE)[uPExDisease_CCLE["Gastrointestinal system disease",]==1], intersect(colnames(uPExMsigDB_CCLE)[uPExMsigDB_CCLE["CHEN_METABOLIC_SYNDROM_NETWORK",]==1], colnames(uPExGOCC_CCLE)[uPExGOCC_CCLE["GO:0031090",]==1]))), intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Diabetes mellitus",]==1], intersect(colnames(uPExDisease_GTEX)[uPExDisease_GTEX["Gastrointestinal system disease",]==1], colnames(uPExMsigDB_GTEX)[uPExMsigDB_GTEX["CHEN_METABOLIC_SYNDROM_NETWORK",]==1]))))

load("../nextProt_pepe.RData") # nextProtXChr
nextprotXensembl <- read.table("../nextprot_ensg_pepe.txt", header=F, sep="\t") # 20350 protes: 17874 PE1 - 1254 uPE1 (1243 with ENSG) / 1899 MPs - 559 
metabolism_uPE_info <- nextProtXChr[nextProtXChr[,2] %in% nextprotXensembl[nextprotXensembl[,2] %in% metabolism_uPE[1],1],]

##########################################
## PE1 DATA
TCGA_uPExPE1 <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_corMatrix.txt", header=T, sep="\t")
rownames(TCGA_uPExPE1) <- TCGA_uPExPE1[,1]
TCGA_uPExPE1 <- TCGA_uPExPE1[,-1]
TCGA_geneXGOBP <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_GeneXenriched_GOBP_filtered.txt", header=T, sep="\t")
rownames(TCGA_geneXGOBP) <- TCGA_geneXGOBP[,1]
TCGA_geneXGOBP <- TCGA_geneXGOBP[,-1]
TCGA_geneXGOCC <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_GeneXenriched_GOCC_filtered.txt", header=T, sep="\t")
rownames(TCGA_geneXGOCC) <- TCGA_geneXGOCC[,1]
TCGA_geneXGOCC <- TCGA_geneXGOCC[,-1]
TCGA_geneXDisease <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_GeneXenriched_Disease_filtered.txt", header=T, sep="\t")
rownames(TCGA_geneXDisease) <- TCGA_geneXDisease[,1]
TCGA_geneXDisease <- TCGA_geneXDisease[,-1]
TCGA_geneXMSigDB  <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_GeneXenriched_MSigDB.txt", header=T, sep="\t")
rownames(TCGA_geneXMSigDB) <- TCGA_geneXMSigDB[,1]
TCGA_geneXMSigDB <- TCGA_geneXMSigDB[,-1]
TCGA_geneXGOMF <- read.table("../../gserranos/PageRank/Data/TCGA/TCGA_GeneXenriched_GOMF_filtered.txt", header=T, sep="\t")
rownames(TCGA_geneXGOMF) <- TCGA_geneXGOMF[,1]
TCGA_geneXGOMF <- TCGA_geneXGOMF[,-1]

CCLE_uPExPE1 <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_corMatrix.txt", header=T, sep="\t")
rownames(CCLE_uPExPE1) <- CCLE_uPExPE1[,1]
CCLE_uPExPE1 <- CCLE_uPExPE1[,-1]
CCLE_geneXGOBP <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_GeneXenriched_GOBP_filtered.txt", header=T, sep="\t")
rownames(CCLE_geneXGOBP) <- CCLE_geneXGOBP[,1]
CCLE_geneXGOBP <- CCLE_geneXGOBP[,-1]
CCLE_geneXGOCC <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_GeneXenriched_GOCC_filtered.txt", header=T, sep="\t")
rownames(CCLE_geneXGOCC) <- CCLE_geneXGOCC[,1]
CCLE_geneXGOCC <- CCLE_geneXGOCC[,-1]
CCLE_geneXDisease <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_GeneXenriched_Disease_filtered.txt", header=T, sep="\t")
rownames(CCLE_geneXDisease) <- CCLE_geneXDisease[,1]
CCLE_geneXDisease <- CCLE_geneXDisease[,-1]
CCLE_geneXMSigDB  <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_GeneXenriched_MSigDB.txt", header=T, sep="\t")
rownames(CCLE_geneXMSigDB) <- CCLE_geneXMSigDB[,1]
CCLE_geneXMSigDB <- CCLE_geneXMSigDB[,-1]
CCLE_geneXGOMF <- read.table("../../gserranos/PageRank/Data/CCLE/CCLE_GeneXenriched_GOMF_filtered.txt", header=T, sep="\t")
rownames(CCLE_geneXGOMF) <- CCLE_geneXGOMF[,1]
CCLE_geneXGOMF <- CCLE_geneXGOMF[,-1]

GTEX_uPExPE1 <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_corMatrix.txt", header=T, sep="\t")
rownames(GTEX_uPExPE1) <- GTEX_uPExPE1[,1]
GTEX_uPExPE1 <- GTEX_uPExPE1[,-1]
GTEX_geneXGOBP <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_GeneXenriched_GOBP_filtered.txt", header=T, sep="\t")
rownames(GTEX_geneXGOBP) <- GTEX_geneXGOBP[,1]
GTEX_geneXGOBP <- GTEX_geneXGOBP[,-1]
GTEX_geneXGOCC <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_GeneXenriched_GOCC_filtered.txt", header=T, sep="\t")
rownames(GTEX_geneXGOCC) <- GTEX_geneXGOCC[,1]
GTEX_geneXGOCC <- GTEX_geneXGOCC[,-1]
GTEX_geneXDisease <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_GeneXenriched_Disease_filtered.txt", header=T, sep="\t")
rownames(GTEX_geneXDisease) <- GTEX_geneXDisease[,1]
GTEX_geneXDisease <- GTEX_geneXDisease[,-1]
GTEX_geneXMSigDB  <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_GeneXenriched_MSigDB.txt", header=T, sep="\t")
rownames(GTEX_geneXMSigDB) <- GTEX_geneXMSigDB[,1]
GTEX_geneXMSigDB <- GTEX_geneXMSigDB[,-1]
GTEX_geneXGOMF <- read.table("../../gserranos/PageRank/Data/GTEX/GTEX_GeneXenriched_GOMF_filtered.txt", header=T, sep="\t")
rownames(GTEX_geneXGOMF) <- GTEX_geneXGOMF[,1]
GTEX_geneXGOMF <- GTEX_geneXGOMF[,-1]

# 1----------- 
## METABOLISM
# TCGA PE genes ENSG00000167483 (720cor-381memb-161diabetes-123Gastro-127metNet) -> 9 PE1
metabolism_TCGA_PE1_1 <- intersect(rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,metabolism_uPE[1]]>0],intersect(rownames(TCGA_geneXMSigDB)[TCGA_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Diabetes.mellitus"]>0],rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0]))))
# TCGA PE genes ENSG00000187808 (924cor) -> 18 PE1
metabolism_TCGA_PE1_2 <- intersect(rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,metabolism_uPE[2]]>0],intersect(rownames(TCGA_geneXMSigDB)[TCGA_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Diabetes.mellitus"]>0],rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0]))))
# TCGA PE genes ENSG00000187862 (1026cor) -> 16 PE1
metabolism_TCGA_PE1_3 <- intersect(rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,metabolism_uPE[3]]>0],intersect(rownames(TCGA_geneXMSigDB)[TCGA_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Diabetes.mellitus"]>0],rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0]))))
# TCGA PE genes ENSG00000204161 (1412cor) -> 27 PE1
metabolism_TCGA_PE1_4 <- intersect(rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,metabolism_uPE[4]]>0],intersect(rownames(TCGA_geneXMSigDB)[TCGA_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Diabetes.mellitus"]>0],rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0]))))
# TCGA PE genes ENSG00000253304 (432cor) -> 4 PE1 - No cruza conlos demás excepto con el 4 (ENSG00000119699)
metabolism_TCGA_PE1_5 <- intersect(rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,metabolism_uPE[5]]>0],intersect(rownames(TCGA_geneXMSigDB)[TCGA_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Diabetes.mellitus"]>0],rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0]))))
metabolism_TCGA_PE1 <- intersect(metabolism_TCGA_PE1_1,intersect(metabolism_TCGA_PE1_2,intersect(metabolism_TCGA_PE1_3,intersect(metabolism_TCGA_PE1_4,metabolism_TCGA_PE1_5))))
# 8 PE1 METABOLISM - TCGA
metabolism_TCGA_PE1 <- intersect(metabolism_TCGA_PE1_1,intersect(metabolism_TCGA_PE1_2,intersect(metabolism_TCGA_PE1_3,metabolism_TCGA_PE1_4)))

# CCLE PE genes ENSG00000167483 (1889cor-1019memb-362diabetes-367Gastro-263metNet) -> 23 PE1
metabolism_CCLE_PE1_1 <- intersect(rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,metabolism_uPE[1]]>0],intersect(rownames(CCLE_geneXMSigDB)[CCLE_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Diabetes.mellitus"]>0],rownames(CCLE_geneXGOCC)[CCLE_geneXGOCC[,"GO.0016020"]>0]))))
# CCLE PE genes ENSG00000187808 (183cor) -> 7 PE1
metabolism_CCLE_PE1_2 <- intersect(rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,metabolism_uPE[2]]>0],intersect(rownames(CCLE_geneXMSigDB)[CCLE_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Diabetes.mellitus"]>0],rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0]))))
# CCLE PE genes ENSG00000187862 (1544cor) -> 16 PE1
metabolism_CCLE_PE1_3 <- intersect(rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,metabolism_uPE[3]]>0],intersect(rownames(CCLE_geneXMSigDB)[CCLE_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Diabetes.mellitus"]>0],rownames(CCLE_geneXGOCC)[CCLE_geneXGOCC[,"GO.0016020"]>0]))))
# CCLE PE genes ENSG00000204161 (694cor) -> 15 PE1
metabolism_CCLE_PE1_4 <- intersect(rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,metabolism_uPE[4]]>0],intersect(rownames(CCLE_geneXMSigDB)[CCLE_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Diabetes.mellitus"]>0],rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0]))))
# CCLE PE genes ENSG00000253304 (2833cor) -> 25 PE1
metabolism_CCLE_PE1_5 <- intersect(rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,metabolism_uPE[5]]>0],intersect(rownames(CCLE_geneXMSigDB)[CCLE_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Diabetes.mellitus"]>0],rownames(CCLE_geneXGOCC)[CCLE_geneXGOCC[,"GO.0016020"]>0]))))
# 5 PE1 METABOLISM - CCLE
metabolism_CCLE_PE1 <- intersect(metabolism_CCLE_PE1_1,intersect(metabolism_CCLE_PE1_2,intersect(metabolism_CCLE_PE1_3,intersect(metabolism_CCLE_PE1_4,metabolism_CCLE_PE1_5))))

# GTEX PE genes ENSG00000167483 (3997cor-1872memb-733diabetes-688Gastro-457metNet) -> 35 PE1
metabolism_GTEX_PE1_1 <- intersect(rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,metabolism_uPE[1]]>0],intersect(rownames(GTEX_geneXMSigDB)[GTEX_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Diabetes.mellitus"]>0],rownames(GTEX_geneXGOCC)[GTEX_geneXGOCC[,"GO.0016020"]>0]))))
# GTEX PE genes ENSG00000187808 (3912cor) -> 32 PE1
metabolism_GTEX_PE1_2 <- intersect(rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,metabolism_uPE[2]]>0],intersect(rownames(GTEX_geneXMSigDB)[GTEX_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Diabetes.mellitus"]>0],rownames(GTEX_geneXGOCC)[GTEX_geneXGOCC[,"GO.0016020"]>0]))))
# GTEX PE genes ENSG00000187862 (1823cor) -> 17 PE1
metabolism_GTEX_PE1_3 <- intersect(rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,metabolism_uPE[3]]>0],intersect(rownames(GTEX_geneXMSigDB)[GTEX_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Diabetes.mellitus"]>0],rownames(GTEX_geneXGOCC)[GTEX_geneXGOCC[,"GO.0016020"]>0]))))
# GTEX PE genes ENSG00000204161 (1895cor) -> 21 PE1
metabolism_GTEX_PE1_4 <- intersect(rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,metabolism_uPE[4]]>0],intersect(rownames(GTEX_geneXMSigDB)[GTEX_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Diabetes.mellitus"]>0],rownames(GTEX_geneXGOCC)[GTEX_geneXGOCC[,"GO.0016020"]>0]))))
# GTEX PE genes ENSG00000253304 (3874cor) -> 35 PE1
metabolism_GTEX_PE1_5 <- intersect(rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,metabolism_uPE[5]]>0],intersect(rownames(GTEX_geneXMSigDB)[GTEX_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Diabetes.mellitus"]>0],rownames(GTEX_geneXGOCC)[GTEX_geneXGOCC[,"GO.0016020"]>0]))))
# 8 PE1 METABOLISM - GTEX
metabolism_GTEX_PE1 <- intersect(metabolism_GTEX_PE1_1,intersect(metabolism_GTEX_PE1_2,intersect(metabolism_GTEX_PE1_3,intersect(metabolism_GTEX_PE1_4,metabolism_GTEX_PE1_5))))

# 4 PE1 in common
intersect(metabolism_TCGA_PE1,intersect(metabolism_CCLE_PE1,metabolism_GTEX_PE1))


# 2----------- 
## IMMUNE
# TCGA PE genes ENSG00000185905 (1249cor-671memb-441autoimmune-115foxp3) -> 32 PE1
immune_TCGA_PE1_1 <- intersect(rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,Immune_uPE[1]]>0],intersect(rownames(TCGA_geneXMSigDB)[TCGA_geneXMSigDB[,"MARSON_BOUND_BY_FOXP3_UNSTIMULATED"]>0],	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Autoimmune.hypersensitivity.disease"]>0]
	,rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0])))
# TCGA PE genes ENSG00000204165 (1259cor) -> 10 PE1
immune_TCGA_PE1_2 <- intersect(rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,Immune_uPE[2]]>0],intersect(rownames(TCGA_geneXMSigDB)[TCGA_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(TCGA_geneXDisease)[TCGA_geneXDisease[,"Diabetes.mellitus"]>0],rownames(TCGA_geneXGOCC)[TCGA_geneXGOCC[,"GO.0016020"]>0]))))
# 0 PE1 IMMUNE - TCGA
immune_TCGA_PE1 <- intersect(immune_TCGA_PE1_1,immune_TCGA_PE1_2)

# CCLE PE genes ENSG00000185905 (1846cor-986memb-285autoimmune-163foxp3) -> 25 PE1
immune_CCLE_PE1_1 <- intersect(rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,Immune_uPE[1]]>0],intersect(rownames(CCLE_geneXMSigDB)[CCLE_geneXMSigDB[,"MARSON_BOUND_BY_FOXP3_UNSTIMULATED"]>0],	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Autoimmune.hypersensitivity.disease"]>0]
	,rownames(CCLE_geneXGOCC)[CCLE_geneXGOCC[,"GO.0016020"]>0])))
# CCLE PE genes ENSG00000204165 (2777cor) -> 27 PE1
immune_CCLE_PE1_2 <- intersect(rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,Immune_uPE[2]]>0],intersect(rownames(CCLE_geneXMSigDB)[CCLE_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(CCLE_geneXDisease)[CCLE_geneXDisease[,"Diabetes.mellitus"]>0],rownames(CCLE_geneXGOCC)[CCLE_geneXGOCC[,"GO.0016020"]>0]))))
# 0 PE1 IMMUNE - CCLE
immune_CCLE_PE1 <- intersect(immune_CCLE_PE1_1,immune_CCLE_PE1_2)

# GTEX PE genes ENSG00000185905 (4683cor-2168memb-701autoimmune-367foxp3) -> 50 PE1
immune_GTEX_PE1_1 <- intersect(rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,Immune_uPE[1]]>0],intersect(rownames(GTEX_geneXMSigDB)[GTEX_geneXMSigDB[,"MARSON_BOUND_BY_FOXP3_UNSTIMULATED"]>0],	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Autoimmune.hypersensitivity.disease"]>0]
	,rownames(GTEX_geneXGOCC)[GTEX_geneXGOCC[,"GO.0016020"]>0])))
# GTEX PE genes ENSG00000204165 (4250cor) -> 27 PE1
immune_GTEX_PE1_2 <- intersect(rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,Immune_uPE[2]]>0],intersect(rownames(GTEX_geneXMSigDB)[GTEX_geneXMSigDB[,"CHEN_METABOLIC_SYNDROM_NETWORK"]>0],	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Gastrointestinal.system.disease"]>0],
	intersect(rownames(GTEX_geneXDisease)[GTEX_geneXDisease[,"Diabetes.mellitus"]>0],rownames(GTEX_geneXGOCC)[GTEX_geneXGOCC[,"GO.0016020"]>0]))))
# 1 PE1 IMMUNE - GTEX
immune_GTEX_PE1 <- intersect(immune_GTEX_PE1_1,immune_GTEX_PE1_2)

# 15 PE1 in common for ENSG00000185905
intersect(immune_TCGA_PE1_1,intersect(immune_CCLE_PE1_1,immune_GTEX_PE1_1))
# 9 PE1 in common for ENSG00000204165
intersect(immune_TCGA_PE1_2,intersect(immune_CCLE_PE1_2,immune_GTEX_PE1_2))

##################################################################################
## Public Experiments (he perdido código: rehago sólo partes necesarias)

library(GEOquery)
library(limma)

genes_sel # Fam129c instead of NIBAN3
genecodev28_Annot
geneXGO
geneXGO_CC
msigdb_mat_ensg2
malaCards_mat_ensg2

#------------------------------------------------
## GSE125746_VAV1.KO.Macrophages_n3_Arrays_Mouse

dataClinic125746
gse125746_cel
gse125746_data
## All samples
gse125746_eset_rma
gse125746_mat_rma
gse125746_mat_rma_filter

design_gse125746
cont.matrix_gse125746
fit2_gse125746

Vav1KOvsWT
write.table(Vav1KOvsWT[Vav1KOvsWT[,2] %in% genes_sel,], file="GSE125746_DEGsel.txt",row.names=F, sep="\t")

sel_gse125746 <- sel
Vav1KOvsWT_enrichedGO
Vav1KOvsWT_enrichedGOcc
Vav1KOvsWT_enrichedMalaCards
Vav1KOvsWT_enrichedMSigDB
cor_tmem273
write.table(gse125746_cor_res, file="gse125746_cor_all.txt", sep="\t", quote=F)

## Vav1 samples
gse125746_data2 <- read.celfiles(gse125746_cel[c(1:3,10:12)])
gse125746_eset_rma2 <- rma(gse125746_data2)
gse125746_mat_rma2 <- exprs(gse125746_eset_rma2)
gse125746_mat_rma_filter2 <- gse125746_mat_rma2[apply(gse125746_mat_rma2[,1:3],1,FUN=function(x) sum(x>5))>=2 | apply(gse125746_mat_rma2[,4:6],1,FUN=function(x) sum(x>5))>=2,]

design_gse125746_2 <- model.matrix(~0+factor(c(rep("Vav1",3),rep("WT",3))))
rownames(design_gse125746_2) <- colnames(gse125746_mat_rma_filter2)
colnames(design_gse125746_2) <- c("Vav1", "WT")
cont.matrix_gse125746_2 <- makeContrasts(Vav1KOvsWT_2 = Vav1 - WT,
	levels=design_gse125746_2)
fit_gse125746 <- lmFit(gse125746_mat_rma_filter2, design_gse125746_2)
fit2_gse125746_2 <- contrasts.fit(fit_gse125746, cont.matrix_gse125746_2)
fit2_gse125746_2 <- eBayes(fit2_gse125746_2)
Vav1KOvsWT_2 <- topTable(fit2_gse125746_2, coef = "Vav1KOvsWT_2", n=nrow(fit2_gse125746_2), adjust="fdr")
Vav1KOvsWT_2 <- cbind(ID=rownames(Vav1KOvsWT_2), Name=sapply(rownames(Vav1KOvsWT_2), FUN=function(x) unlist(geneName125746[x])), Vav1KOvsWT_2)
Vav1KOvsWT_2[Vav1KOvsWT_2[,2] %in% genes_sel,]
write.table(Vav1KOvsWT_2[Vav1KOvsWT_2[,2] %in% genes_sel,], file="GSE125746_DEGsel_Vav1.txt",row.names=F, sep="\t")

#------------------------------------------------
## GSE102889_WAS.KO.TcellLymphoma_n5_Arrays_Mouse

gse102889 <- getGEO(filename="GSE102889/GSE102889_family.soft.gz")
charSamples102889 <- list(length(GSMList(gse102889)))
charDescription102889 <- list(length(GSMList(gse102889)))
charTitle102889 <- list(length(GSMList(gse102889)))
for (i in 1:length(GSMList(gse102889))) {                                
	charSamples102889[i] <- list(Meta(GSMList(gse102889)[[i]])$characteristics_ch1)
}
for (i in 1:length(GSMList(gse102889))) {                                
	charDescription102889[i] <- list(Meta(GSMList(gse102889)[[i]])$geo_accession)
}
for (i in 1:length(GSMList(gse102889))) {                                
	charTitle102889[i] <- list(Meta(GSMList(gse102889)[[i]])$title)
}
dataClinic102889 <- data.frame(GEO=unlist(charDescription102889),Sample=unlist(charTitle102889), Tissue= sapply(strsplit(sapply(charSamples102889,FUN=function(x) x[1]),": "), FUN=function(x) x[2]), Genotype= sapply(strsplit(sapply(charSamples102889,FUN=function(x) x[2]),": "), FUN=function(x) x[2]))
rownames(dataClinic102889) <- dataClinic102889[,1]
load("~/dato-activo/00_References/AffyAnnot/MoGene10stv1_na35_mm10.rda")

gse102889_cel <- list.files(path="GSE102889",pattern=".CEL.gz", full.names=T)
gse102889_data <- read.celfiles(gse102889_cel)
gse102889_eset_rma <- rma(gse102889_data)
gse102889_mat_rma <- exprs(gse102889_eset_rma)

tipo102889 <- factor(substr(paste(dataClinic102889[substr(colnames(gse102889_mat_rma),1,10),"Genotype"]),1,2))
gse102889_mat_rma_filter <- gse102889_mat_rma[apply(gse102889_mat_rma[,tipo102889=="WT"],1,FUN=function(x) sum(x>5))>=3 | apply(gse102889_mat_rma[,tipo102889=="WA"],1,FUN=function(x) sum(x>5))>=3,]
load("~/dato-activo/00_References/AffyAnnot/HuGene20stv1_na36_hg19.rda")
gse83452_mat_rma_filter_annot <- merge(annot2,data.frame(Probeset=rownames(gse102889_mat_rma_filter), gse102889_mat_rma_filter), by.x = 1, by.y = 1)

# QC
acol = sample(brewer.pal(8, "Dark2"), ncol(gse102889_mat_rma), replace = (8 < ncol(gse102889_mat_rma)))
xlim_r = c(min(na.omit(exprs(gse102889_data))),max(na.omit(exprs(gse102889_data)))) 
xlim = c(min(na.omit(gse102889_mat_rma)),max(na.omit(gse102889_mat_rma)))
clust.euclid.average <- hclust(dist(t(exprs(gse102889_data))),method="average")
clust.euclid.average2 <- hclust(dist(t(gse102889_mat_rma)),method="average")
pdf(file="GSE102889_QC_Array_Mar20.pdf", width = 10, height = 10)
boxplot(exprs(gse102889_data),names=colnames(gse102889_eset_rma), which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", exprs(gse102889_data), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(gse102889_mat_rma, cex.axis=0.7, las=2)
multi("density", gse102889_mat_rma, xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()
a <- detOutliers(gse102889_mat_rma, tipo102889, "GSE102889_Array_var_Mar20.pdf", 20, 20)

x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(gse102889_mat_rma))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo102889)) + geom_text(aes(label=colnames(gse102889_mat_rma)), hjust=0, vjust=0, nudge_x=-8, nudge_y=2) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "GSE102889_PCA_Mar20.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

design_gse102889 <- model.matrix(~0+tipo102889)
rownames(design_gse102889) <- colnames(gse102889_mat_rma_filter)
colnames(design_gse102889) <- levels(tipo102889)
cont.matrix_gse102889 <- makeContrasts(WASPvsWT = WA - WT,
	levels=design_gse102889)
fit_gse102889 <- lmFit(gse102889_mat_rma_filter, design_gse102889)
fit2_gse102889 <- contrasts.fit(fit_gse102889, cont.matrix_gse102889)
fit2_gse102889 <- eBayes(fit2_gse102889)
WASPvsWT <- topTable(fit2_gse102889, coef = "WASPvsWT", n=nrow(fit2_gse102889), adjust="fdr")
WASPvsWT <- merge(annot2,data.frame(Probeset=rownames(WASPvsWT), WASPvsWT), by.x = 1, by.y = 1)
WASPvsWT <- WASPvsWT[order(WASPvsWT$P.Value),]
WASPvsWT[WASPvsWT[,3] %in% genes_sel,]
write.table(WASPvsWT[WASPvsWT[,3] %in% genes_sel,], file="GSE102889_DEGsel_Was.txt",row.names=F, sep="\t")

cor <- apply(gse102889_mat_rma_filter[paste(WASPvsWT[WASPvsWT[,3] %in% genes_sel,1]),], 1, FUN=function(x) cor(as.numeric(paste(gse102889_mat_rma_filter["10508707",])),as.numeric(paste(x))))
p <- apply(gse102889_mat_rma_filter[paste(WASPvsWT[WASPvsWT[,3] %in% genes_sel,1]),], 1, FUN=function(x) cor.test(as.numeric(paste(gse102889_mat_rma_filter["10508707",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")
gse102889followUp_cor_res <- data.frame(ID=names(cor),gene=WASPvsWT[names(cor),3],cor=cor, p=p, fdr=fdr)
write.table(gse102889followUp_cor_res, file="GSE102889_Tmem200b_cor.txt", row.names=F, sep="\t", quote=F)

#------------------------------------------------
## GSE109022_InsulinTreatment_n4n3_Arrays_Human

gse109022 <- getGEO(filename="GSE109022/GSE109022_family.soft.gz")
charSamples109022 <- list(length(GSMList(gse109022)))
charDescription109022 <- list(length(GSMList(gse109022)))
charTitle109022 <- list(length(GSMList(gse109022)))
for (i in 1:length(GSMList(gse109022))) {                                
	charSamples109022[i] <- list(Meta(GSMList(gse109022)[[i]])$characteristics_ch1)
}
for (i in 1:length(GSMList(gse109022))) {                                
	charDescription109022[i] <- list(Meta(GSMList(gse109022)[[i]])$geo_accession)
}
for (i in 1:length(GSMList(gse109022))) {                                
	charTitle109022[i] <- list(Meta(GSMList(gse109022)[[i]])$title)
}
dataClinic109022 <- data.frame(GEO=unlist(charDescription109022),Sample=unlist(charTitle109022), Treatment= sapply(strsplit(sapply(charSamples109022,FUN=function(x) x[1]),": "), FUN=function(x) x[2]), CellLine= sapply(strsplit(sapply(charSamples109022,FUN=function(x) x[2]),": "), FUN=function(x) x[2]))
rownames(dataClinic109022) <- dataClinic109022[,1]
annotHuGene <- read.table(file = "GSE109022/HuGene-2_1-st-v1.na36.hg19.transcript.csv", header=T, sep = ",")
annot <- data.frame("probeset_id"=annotHuGene$transcript_cluster_id, parseAnnotGene(annotHuGene[,8:9]))
ann <- annot[-grep("_control", annotHuGene[,9]),]
save(list="ann",file="~/dato-activo/00_References/AffyAnnot/HuGene21stv1_na36_hg19.rda")

gse109022_cel <- list.files(path="GSE109022",pattern=".cel.gz", full.names=T)
gse109022_data <- read.celfiles(gse109022_cel)
gse109022_eset_rma <- rma(gse109022_data)
gse109022_mat_rma <- exprs(gse109022_eset_rma)
colnames(gse109022_mat_rma) <- gsub(".ga.cel.gz","",colnames(gse109022_mat_rma))

tipo109022 <- substr(paste(dataClinic109022[substr(colnames(gse109022_mat_rma),1,10),"Treatment"]),14,15)
tipo109022[-grep("[d-i]",tipo109022)] <- "C"
tipo109022[grep("d",tipo109022)] <- "Detemir"
tipo109022[grep("-l",paste(dataClinic109022[substr(colnames(gse109022_mat_rma),1,10),"Treatment"]))] <- "IGF1"
tipo109022[grep("gla",paste(dataClinic109022[substr(colnames(gse109022_mat_rma),1,10),"Treatment"]))] <- "Glargine"
tipo109022[grep("in$",tipo109022)] <- "Insulin"
tipo109022 <- factor(tipo109022)
cell109022 <- factor(gsub("-","",paste(dataClinic109022[substr(colnames(gse109022_mat_rma),1,10),"CellLine"])))
gse109022_mat_rma_filter <- gse109022_mat_rma[apply(gse109022_mat_rma[,tipo109022=="C" & cell109022=="ECC1"],1,FUN=function(x) sum(x>4))>=2 | apply(gse109022_mat_rma[,tipo109022=="IGF1" & cell109022=="ECC1"],1,FUN=function(x) sum(x>4))>=2 | apply(gse109022_mat_rma[,tipo109022=="Insulin" & cell109022=="ECC1"],1,FUN=function(x) sum(x>4))>=2,]
gse109022_mat_rma_filter_annot <- data.frame(ann[rownames(gse109022_mat_rma_filter),],gse109022_mat_rma_filter)

# QC
acol = sample(brewer.pal(8, "Dark2"), ncol(gse109022_mat_rma), replace = (8 < ncol(gse109022_mat_rma)))
xlim_r = c(min(na.omit(log2(exprs(gse109022_data)+1))),max(na.omit(log2(exprs(gse102889_data)))) )
xlim = c(min(na.omit(gse109022_mat_rma)),max(na.omit(gse109022_mat_rma)))
clust.euclid.average <- hclust(dist(t(log2(exprs(gse109022_data)+1))),method="average")
clust.euclid.average2 <- hclust(dist(t(gse109022_mat_rma)),method="average")
pdf(file="GSE109022_QC_Array_Mar20.pdf", width = 10, height = 10)
boxplot(log2(exprs(gse109022_data)+1),names=colnames(gse109022_eset_rma), which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(exprs(gse109022_data)), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(gse109022_mat_rma, cex.axis=0.7, las=2)
multi("density", gse109022_mat_rma, xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()
a <- detOutliers(gse109022_mat_rma, tipo109022, "GSE109022_Array_varTipo_Abr20.pdf", 10, 10)
a <- detOutliers(gse109022_mat_rma, cell109022, "GSE109022_Array_varCell_Abr20.pdf", 10, 10)

x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(gse109022_mat_rma))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo109022, shape=cell109022)) + geom_text(aes(label=colnames(gse109022_mat_rma)), hjust=0, vjust=0, nudge_x=-8, nudge_y=2) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "GSE109022_PCA_Abr2020.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

design_gse109022 <- model.matrix(~0+factor(paste(tipo109022,cell109022,sep="_")))
rownames(design_gse109022) <- colnames(gse109022_mat_rma_filter)
colnames(design_gse109022) <- levels(factor(paste(tipo109022,cell109022,sep="_")))
cont.matrix_gse109022 <- makeContrasts(InsulinVsC = Insulin_ECC1 - C_ECC1,
	IGF1vsC = IGF1_ECC1 - C_ECC1,
	levels=design_gse109022)
fit_gse109022 <- lmFit(gse109022_mat_rma_filter, design_gse109022)
fit2_gse109022 <- contrasts.fit(fit_gse109022, cont.matrix_gse109022)
fit2_gse109022 <- eBayes(fit2_gse109022)
InsulinVsC <- topTable(fit2_gse109022, coef = "InsulinVsC", n=nrow(fit2_gse109022), adjust="fdr")
InsulinVsC <- data.frame(ann[rownames(InsulinVsC),],InsulinVsC)
IGF1vsC <- topTable(fit2_gse109022, coef = "IGF1vsC", n=nrow(fit2_gse109022), adjust="fdr")
IGF1vsC <- data.frame(ann[rownames(IGF1vsC),],IGF1vsC)
write.table(InsulinVsC[InsulinVsC[,3] %in% toupper(genes_sel),], file="GSE109022_DEGsel_InsulinVsC.txt",row.names=F, sep="\t")
write.table(IGF1vsC[IGF1vsC[,3] %in% toupper(genes_sel),], file="GSE109022_DEGsel_IGF1vsC.txt",row.names=F, sep="\t")

#------------------------------------------------
## GSE95489_MetabolicSyndromeModel_n3_Arrays_Mouse

gse95489 <- getGEO(filename="GSE95489/GSE95489_family.soft.gz")
charSamples95489 <- list(length(GSMList(gse95489)))
charDescription95489 <- list(length(GSMList(gse95489)))
charTitle95489 <- list(length(GSMList(gse95489)))
for (i in 1:length(GSMList(gse95489))) {                                
	charSamples95489[i] <- list(Meta(GSMList(gse95489)[[i]])$characteristics_ch1)
}
for (i in 1:length(GSMList(gse95489))) {                                
	charDescription95489[i] <- list(Meta(GSMList(gse95489)[[i]])$geo_accession)
}
for (i in 1:length(GSMList(gse95489))) {                                
	charTitle95489[i] <- list(Meta(GSMList(gse95489)[[i]])$title)
}
dataClinic95489 <- data.frame(GEO=unlist(charDescription95489),Sample=unlist(charTitle95489), Gender= sapply(strsplit(sapply(charSamples95489,FUN=function(x) x[1]),": "), FUN=function(x) x[2]), Age= sapply(strsplit(sapply(charSamples95489,FUN=function(x) x[2]),": "), FUN=function(x) x[2]), Strain= sapply(strsplit(sapply(charSamples95489,FUN=function(x) x[3]),": "), FUN=function(x) x[2]), Tissue= sapply(strsplit(sapply(charSamples95489,FUN=function(x) x[4]),": "), FUN=function(x) x[2]))
rownames(dataClinic95489) <- dataClinic95489[,1]
load("~/dato-activo/00_References/AffyAnnot/MoGene10stv1_na35_mm10.rda")

gse95489_cel <- list.files(path="GSE95489",pattern=".CEL.gz", full.names=T)
gse95489_data <- read.celfiles(gse95489_cel)
gse95489_eset_rma <- rma(gse95489_data)
gse95489_mat_rma <- exprs(gse95489_eset_rma)
colnames(gse95489_mat_rma) <- gsub(".CEL.gz","",colnames(gse95489_mat_rma))

tipo95489 <- factor(c(rep("WT",3),rep("NZ0",3)))
gse95489_mat_rma_filter <- gse95489_mat_rma[apply(gse95489_mat_rma[,tipo95489=="WT"],1,FUN=function(x) sum(x>4))>=2 | apply(gse95489_mat_rma[,tipo95489=="NZ0"],1,FUN=function(x) sum(x>4))>=2,]
library(mouse4302.db)
gse95489_mat_rma_filter_annot <- data.frame(Probe=rownames(gse95489_mat_rma_filter),GeneName=mapIds(mouse4302.db,keys=rownames(gse95489_mat_rma_filter),keytype="PROBEID",column="SYMBOL",multiVals="first"),gse95489_mat_rma_filter)

# QC
acol = sample(brewer.pal(8, "Dark2"), ncol(gse95489_mat_rma), replace = (8 < ncol(gse95489_mat_rma)))
xlim_r = c(min(na.omit(log2(exprs(gse95489_data)))),max(na.omit(log2(exprs(gse95489_data))))) 
xlim = c(min(na.omit(gse95489_mat_rma)),max(na.omit(gse95489_mat_rma)))
clust.euclid.average <- hclust(dist(t(log2(exprs(gse95489_data)))),method="average")
clust.euclid.average2 <- hclust(dist(t(gse95489_mat_rma)),method="average")
pdf(file="GSE95489_QC_Array_Abr20.pdf", width = 10, height = 10)
boxplot(log2(exprs(gse95489_data)),names=colnames(gse95489_eset_rma), which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(exprs(gse95489_data)), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(gse95489_mat_rma, cex.axis=0.7, las=2)
multi("density", gse95489_mat_rma, xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()
a <- detOutliers(gse95489_mat_rma, tipo95489, "GSE95489_Array_var_Abr20.pdf", 10, 10)

x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(gse95489_mat_rma))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo95489)) + geom_text(aes(label=colnames(gse95489_mat_rma)), hjust=0, vjust=0, nudge_x=-8, nudge_y=2) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "GSE95489_PCA_Mar20.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

design_gse95489 <- model.matrix(~0+tipo95489)
rownames(design_gse95489) <- colnames(gse95489_mat_rma_filter_annot[,-c(1:2)])
colnames(design_gse95489) <- levels(tipo95489)
cont.matrix_gse95489 <- makeContrasts(NZ0vsWT = NZ0 - WT,
	levels=design_gse95489)
fit_gse95489 <- lmFit(gse95489_mat_rma_filter_annot[,-c(1:2)], design_gse95489)
fit2_gse95489 <- contrasts.fit(fit_gse95489, cont.matrix_gse95489)
fit2_gse95489 <- eBayes(fit2_gse95489)
NZ0vsWT <- topTable(fit2_gse95489, coef = "NZ0vsWT", n=nrow(fit2_gse95489), adjust="fdr")
NZ0vsWT <- data.frame(Probe=rownames(NZ0vsWT),GeneName=mapIds(mouse4302.db,keys=rownames(NZ0vsWT),keytype="PROBEID",column="SYMBOL",multiVals="first"),NZ0vsWT)
NZ0vsWT[NZ0vsWT[,2] %in% genes_sel,]
write.table(NZ0vsWT[NZ0vsWT[,2] %in% genes_sel,], file="GSE95489_DEGsel_NZ0vsWT.txt",row.names=F, sep="\t")

cor <- apply(gse95489_mat_rma_filter[paste(NZ0vsWT[NZ0vsWT[,2] %in% genes_sel,1]),], 1, FUN=function(x) cor(as.numeric(paste(gse95489_mat_rma_filter["1429604_at",])),as.numeric(paste(x))))
p <- apply(gse95489_mat_rma_filter[paste(NZ0vsWT[NZ0vsWT[,2] %in% genes_sel,1]),], 1, FUN=function(x) cor.test(as.numeric(paste(gse95489_mat_rma_filter["1429604_at",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")
gse95489_cor_res <- data.frame(ID=names(cor),gene=NZ0vsWT[names(cor),2],cor=cor, p=p, fdr=fdr)
write.table(gse95489_cor_res, file="GSE95489_Tmem273_cor.txt", row.names=F, sep="\t", quote=F)

#----------------------------------------------------
## GSE104674_DiabetesAdiposeVsNormal_n6_RNASeq_Human

dataClinic104674
GSE104674_raw
GSE104674_raw_filter

## All samples
GSE104674_f_norm

design_gse104674
cont.matrix_gse104674
fit2_gse104674

D15vsH15
D545vsH545
D1115vsH1115
D1530vsH1530
gse104674_bmat <- cbind(D15vsH15[rownames(D15vsH15)[D15vsH15[,1] %in% toupper(genes_sel)],c(1:2,5,6)],D545vsH545[rownames(D15vsH15)[D15vsH15[,1] %in% toupper(genes_sel)],c(2,5,6)],D1115vsH1115[rownames(D15vsH15)[D15vsH15[,1] %in% toupper(genes_sel)],c(2,5,6)],D1530vsH1530[rownames(D15vsH15)[D15vsH15[,1] %in% toupper(genes_sel)],c(2,5,6)])
colnames(gse104674_bmat)[-1] <- paste(c(rep("D15",3),rep("D545",3),rep("D1115",3),rep("D1530",3)),colnames(gse104674_bmat)[-1],sep=".")
write.table(gse104674_bmat, file="GSE104674_DEGsel_all.txt",row.names=F, sep="\t")

cor_Fam129c
cor_tmem200b
gse104674_cor_res1
gse104674_cor_res2
write.table(gse104674_bmat, file="GSE104674_cor_all.txt",row.names=F, sep="\t")

## D1530 samples
GSE104674_raw_D1530 <- GSE104674_raw[whichFilter.D1530 | whichFilter.H1530, tipo23 %in% c("Diabetes_1530","Healthy_1530")]
design_gse104674_D1530 <- model.matrix(~0+factor(paste(tipo23[tipo23 %in% c("Diabetes_1530","Healthy_1530")])))
colnames(design_gse104674_D1530) <- c("Diabetes_1530","Healthy_1530")
rownames(design_gse104674_D1530) <- colnames(GSE104674_raw_D1530)
GSE104674_D1530_fn <- calcNormFactors(DGEList(counts=GSE104674_raw_D1530))
GSE104674_D1530_norm <- voom(GSE104674_D1530_fn, design_gse104674_D1530)
cont.matrix_gse104674_D1530 <- makeContrasts(D1530vsH1530_2 = Diabetes_1530 - Healthy_1530,
		levels=design_gse104674_D1530)
fit_gse104674_1530 <- lmFit(GSE104674_D1530_norm, design_gse104674_D1530)
fit2_gse104674_1530 <- contrasts.fit(fit_gse104674_1530, cont.matrix_gse104674_D1530)
fit2_gse104674_1530 <- eBayes(fit2_gse104674_1530)
D1530vsH1530_2 <- topTable(fit2_gse104674_1530, coef = "D1530vsH1530_2", n=nrow(fit2_gse104674_1530), adjust="fdr")
D1530vsH1530_2 <- data.frame(Name=sapply(rownames(D1530vsH1530_2), FUN=function(x) unlist(strsplit(x,":"))[1]), D1530vsH1530_2)
D1530vsH1530_2[D1530vsH1530_2[,1] %in% toupper(genes_sel),]
write.table(D1530vsH1530_2[D1530vsH1530_2[,1] %in% toupper(genes_sel),], file="GSE104674_DEGsel_D1530.txt",row.names=F, sep="\t")

cor <- apply(GSE104674_D1530_norm$E[rownames(D1530vsH1530_2)[D1530vsH1530_2[,1] %in% toupper(genes_sel)],], 1, FUN=function(x) cor(as.numeric(paste(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",])),as.numeric(paste(x))))
p <- apply(GSE104674_D1530_norm$E[rownames(D1530vsH1530_2)[D1530vsH1530_2[,1] %in% toupper(genes_sel)],], 1, FUN=function(x) cor.test(as.numeric(paste(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")
gse104674_cor_res2_D1530 <- data.frame(ID=names(cor),gene=D1530vsH1530_2[names(cor),1],cor=cor, p=p, fdr=fdr)
write.table(gse104674_cor_res, file="GSE104674_cor_allSamples.txt.txt", row.names=F, sep="\t", quote=F)
write.table(gse104674_cor_res2_D1530, file="GSE104674_cor_D1530.txt.txt", row.names=F, sep="\t", quote=F)

pdf("GSE104674_D1530_FAM129C_cor_Mar20.pdf", width = 12, height = 8)
{ 
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",]),"FAM129C_log2exp_2"=as.numeric(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",]))
	p1=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=FAM129C_log2exp_2)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("FAM129C: cor=", round(gse104674_cor_res2_D1530$cor[which(gse104674_cor_res2_D1530[,2]== "FAM129C")],3)," - p=", round(gse104674_cor_res2_D1530$p[which(gse104674_cor_res2_D1530[,2] == "FAM129C")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"FAM129C_log2exp_2"]), hjust=2, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",]),"TMEM200B_log2exp"=as.numeric(GSE104674_D1530_norm$E["TMEM200B:chr1:NM_001171868,NM_001003682",]))
	p2=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=TMEM200B_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("TMEM200B: cor=", round(gse104674_cor_res2_D1530$cor[which(gse104674_cor_res2_D1530[,2]== "TMEM200B")],3)," - p=", round(gse104674_cor_res2_D1530$p[which(gse104674_cor_res2_D1530[,2] == "TMEM200B")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"TMEM200B_log2exp"]), hjust=1.2, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",]),"WAS_log2exp"=as.numeric(GSE104674_D1530_norm$E["WAS:chrX:NM_000377",]))
	p3=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=WAS_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("WAS: cor=", round(gse104674_cor_res2_D1530$cor[which(gse104674_cor_res2_D1530[,2]== "WAS")],3)," - p=", round(gse104674_cor_res2_D1530$p[which(gse104674_cor_res2_D1530[,2] == "WAS")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"WAS_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",]),"VAV1_log2exp"=as.numeric(GSE104674_D1530_norm$E["VAV1:chr19:NM_005428",]))
	p4=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=VAV1_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("VAV1: cor=", round(gse104674_cor_res2_D1530$cor[which(gse104674_cor_res2_D1530[,2]== "VAV1")],3)," - p=", round(gse104674_cor_res2_D1530$p[which(gse104674_cor_res2_D1530[,2] == "VAV1")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"VAV1_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",]),"PLEK_log2exp"=as.numeric(GSE104674_D1530_norm$E["PLEK:chr2:NM_002664",]))
	p5=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=PLEK_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("PLEK: cor=", round(gse104674_cor_res2_D1530$cor[which(gse104674_cor_res2_D1530[,2]== "PLEK")],3)," - p=", round(gse104674_cor_res2_D1530$p[which(gse104674_cor_res2_D1530[,2] == "PLEK")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"PLEK_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(GSE104674_D1530_norm$E["FAM129C:chr19:NM_001098524,NM_173544",]),"INPP5D_log2exp"=as.numeric(GSE104674_D1530_norm$E["INPP5D:chr2:NM_001017915,NM_005541",]))
	p6=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=INPP5D_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("INPP5D: cor=", round(gse104674_cor_res2_D1530$cor[which(gse104674_cor_res2_D1530[,2]== "INPP5D")],3)," - p=", round(gse104674_cor_res2_D1530$p[which(gse104674_cor_res2_D1530[,2] == "INPP5D")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"INPP5D_log2exp"]), hjust=1.4, size=4, vjust=0) + theme_bw()
	multiplot(p1,p2,p3,p4,p5,p6,cols=3)
}
dev.off()
# busco gene-set de genes correlados para FAM129
FAM129C_corTCGA <- rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,substr(rownames(genecodev28_Annot)[genecodev28_Annot$gene_name=="FAM129C"],1,15)]>0]
FAM129C_corCCLE <- rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,substr(rownames(genecodev28_Annot)[genecodev28_Annot$gene_name=="FAM129C"],1,15),]>0]
FAM129C_corGTEX <- rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,substr(rownames(genecodev28_Annot)[genecodev28_Annot$gene_name=="FAM129C"],1,15),]>0]
FAM129C_cor <- intersect(FAM129C_corTCGA,intersect(FAM129C_corCCLE,FAM129C_corGTEX))
GSE104674_FAM129C_phyper <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(D1530vsH1530[,1])), paste(D1530vsH1530[D1530vsH1530$P.Value<0.05,1]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(D1530vsH1530[,1]))), length(unique(D1530vsH1530[,1]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(D1530vsH1530[,1]))), length(unique(D1530vsH1530[D1530vsH1530$P.Value<0.05,1])), lower.tail=F) # 0.99
GSE104674_FAM129C_phyper2 <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(D1530vsH1530[,1])), paste(D1530vsH1530[D1530vsH1530$P.Value<0.01,1]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(D1530vsH1530[,1]))), length(unique(D1530vsH1530[,1]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(D1530vsH1530[,1]))), length(unique(D1530vsH1530[D1530vsH1530$P.Value<0.01,1])), lower.tail=F) # 0.85
pdf("GSE104674_FAM129Cenrichment.pdf",16,8)
par(mfrow=c(1,2))
compare2List(paste(unique(D1530vsH1530[D1530vsH1530$P.Value<0.05,1])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(D1530vsH1530[,1])), "Diabetes1530vsHealthy1530_p05","FAM129Ccor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE104674_FAM129C_phyper,3),sep=""))
compare2List(paste(unique(D1530vsH1530[D1530vsH1530$P.Value<0.01,1])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(D1530vsH1530[,1])), "Diabetes1530vsHealthy1530_p01","FAM129Ccor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE104674_FAM129C_phyper2,3),sep=""))
dev.off()

## D015 samples
GSE104674_raw_D015 <- GSE104674_raw[whichFilter.D015 | whichFilter.H015, tipo23 %in% c("Diabetes_015","Healthy_015")]
design_gse104674_D015 <- model.matrix(~0+factor(paste(tipo23[tipo23 %in% c("Diabetes_015","Healthy_015")])))
colnames(design_gse104674_D015) <- c("Diabetes_015","Healthy_015")
rownames(design_gse104674_D015) <- colnames(GSE104674_raw_D015)
GSE104674_D015_fn <- calcNormFactors(DGEList(counts=GSE104674_raw_D015))
GSE104674_D015_norm <- voom(GSE104674_D015_fn, design_gse104674_D015)
cont.matrix_gse104674_D015 <- makeContrasts(D15vsH15_2 = Diabetes_015 - Healthy_015,
		levels=design_gse104674_D015)
fit_gse104674_015 <- lmFit(GSE104674_D015_norm, design_gse104674_D015)
fit2_gse104674_015 <- contrasts.fit(fit_gse104674_015, cont.matrix_gse104674_D015)
fit2_gse104674_015 <- eBayes(fit2_gse104674_015)
D15vsH15_2 <- topTable(fit2_gse104674_015, coef = "D15vsH15_2", n=nrow(fit2_gse104674_1530), adjust="fdr")
D15vsH15_2 <- data.frame(Name=sapply(rownames(D15vsH15_2), FUN=function(x) unlist(strsplit(x,":"))[1]), D15vsH15_2)
D15vsH15_2[D15vsH15_2[,1] %in% toupper(genes_sel),]


#---------------------------------------------------
## GSE139929_InsulinResistanceVsNormalHepG2_n3_RNASeq_Human
# sbatch --depend=aftercorr:123 my.job

gse139929 <- getGEO(filename="GSE139929_family.soft.gz")
charSamples139929 <- list(length(GSMList(gse139929)))
charDescription139929 <- list(length(GSMList(gse139929)))
charTitle139929 <- list(length(GSMList(gse139929)))
for (i in 1:length(GSMList(gse139929))) {                                
	charSamples139929[i] <- list(Meta(GSMList(gse139929)[[i]])$characteristics_ch1)
}
for (i in 1:length(GSMList(gse139929))) {                                
	charDescription139929[i] <- list(Meta(GSMList(gse139929)[[i]])$geo_accession)
}
for (i in 1:length(GSMList(gse139929))) {                                
	charTitle139929[i] <- list(Meta(GSMList(gse139929)[[i]])$title)
}
dataClinic139929 <- data.frame(GEO=unlist(charDescription139929),Sample=unlist(charTitle139929), Tissue= sapply(strsplit(sapply(charSamples139929,FUN=function(x) x[1]),": "), FUN=function(x) x[2]))
rownames(dataClinic139929) <- dataClinic139929[,1]

genecodev32_all <- read.table(file = "~/dato-activo/00_References/gencode_gtf/gencode.v32.annotation.gtf", skip = 5, header = FALSE, sep = "\t")
genecodev32 <- genecodev32_all[genecodev32_all$V3 == "gene", ]
colnames(genecodev32) <- c("chr","DB","Type", "start","end","","strand","","Description")
genecodev32_tmp <- apply(as.data.frame(genecodev32[,9]), 1, parseENCODE_eli)   
genecodev32_tmp <- t(genecodev32_tmp)
colnames(genecodev32_tmp) <- c("gene_id", "gene_type", "gene_name", "level", "havana_gene")
genecodev32_Annot <- data.frame(genecodev32_tmp, genecodev32[,-c(6,8:9)])
rownames(genecodev32_Annot) <- genecodev32_Annot[,1]
geneDesc_Ensembl <- read.table(file = "~/dato-activo/00_References/gencode_gtf/EnsemblGeneDescription_H98.txt", header = TRUE, sep = "\t")
genecodev32_Annot <- merge(geneDesc_Ensembl, genecodev32_Annot, by.x=2, by.y=1, all.x=F, all.y=T)
rownames(genecodev32_Annot) <- genecodev32_Annot[,1]
colnames(genecodev32_Annot)[1] <- "gene_id"
genecodev32_Annot <- genecodev32_Annot[,c(1,3:4,2,5:ncol(genecodev32_Annot))]

countsDir <- list.files(path="FASTQ",full.names=T)
countsFiles <- list.files(pattern=".counts.txt$", path=countsDir[1], full.names=T)
for(i in 2:length(countsDir))
{
	countsFiles <- c(countsFiles,list.files(pattern=".counts.txt$", path=countsDir[i], full.names=T))
}
dataMatCounts <- data.frame()
for (i in 1:length(countsFiles)) {		

	cat(i, "\n")

	tmp <- read.csv(file = countsFiles[i], header = TRUE, skip = 1, sep = "\t", fill = TRUE)
	colnames(tmp)[7] <- countsFiles[i]

	if (i == 1) {

		dataMatCounts <- tmp

	} else {

		dataMatCounts <- cbind(dataMatCounts, tmp[7])

	}

}
colnames(dataMatCounts) <- unlist(strsplit(gsub("_Clean_Data.counts.txt", "", colnames(dataMatCounts)),"/"))[c(1:6, seq(from=9, to=((ncol(dataMatCounts)-6)*3)+6,by=3))]
rownames(dataMatCounts) <- dataMatCounts[,1]
gse139929_mat <- dataMatCounts
tipo139929 <- factor(c(rep("CTRL",3), rep("MODEL",3)))

# FILTERING 
gse139929_mat_f <- gse139929_mat[apply(gse139929_mat[,-c(1:6)][,tipo139929=="CTRL"], 1, FUN=function(x) sum(x>5))>1 | apply(gse139929_mat[,-c(1:6)][,tipo139929=="MODEL"], 1, FUN=function(x) sum(x>5))>1,]
# Normalization
design139929 <- model.matrix(~0+tipo139929)
colnames(design139929) <- levels(tipo139929)
rownames(design139929) <- colnames(gse139929_mat_f)[-c(1:6)]
dataMatCounts_fn139929 <- calcNormFactors(DGEList(counts=gse139929_mat_f[,-c(1:6)]))
dataMatCounts_norm139929 <- voom(dataMatCounts_fn139929, design139929)
dataMatCounts_norm139929_annot <- data.frame(genecodev32_Annot[rownames(dataMatCounts_norm139929),], dataMatCounts_norm139929)

#######################
# QC
acol = sample(brewer.pal(8, "Dark2"), ncol(gse139929_mat_f)-6, replace = (8 < ncol(gse139929_mat_f)-6))
xlim_r = c(min(na.omit(log2(gse139929_mat_f[,-c(1:6)]+1))),max(na.omit(log2(gse139929_mat_f[,-c(1:6)]+1)))) 
xlim = c(min(na.omit(dataMatCounts_norm139929$E)),max(na.omit(dataMatCounts_norm139929$E)))
clust.euclid.average <- hclust(dist(t(log2(gse139929_mat_f[,-c(1:6)]+1))),method="average")
clust.euclid.average2 <- hclust(dist(t(dataMatCounts_norm139929$E)),method="average")
pdf(file="GSE139929_QC_RNASeq_Mar20.pdf", width = 10, height = 10)
boxplot(log2(gse139929_mat_f[,-c(1:6)]+1),names=colnames(gse139929_mat_f)[-c(1:6)], which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", log2(gse139929_mat_f[,-c(1:6)]+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(dataMatCounts_norm139929$E, cex.axis=0.7, las=2)
multi("density", dataMatCounts_norm139929$E, xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()
a <- detOutliers(dataMatCounts_norm139929$E, tipo139929, "GSE139929_RNASeq_var_Mar20.pdf", 10, 10)

x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(dataMatCounts_norm139929$E))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo139929)) + geom_text(aes(label=colnames(dataMatCounts_norm139929$E)), hjust=0, vjust=0, nudge_x=-8, nudge_y=2) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "GSE139929_PCA_Mar20.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

cont.matrix_gse139929 <- makeContrasts(MODELvsCTRL = MODEL - CTRL,
		levels=design139929)
fit_gse139929 <- lmFit(dataMatCounts_norm139929, design139929)
fit2_gse139929 <- contrasts.fit(fit_gse139929, cont.matrix_gse139929)
fit2_gse139929 <- eBayes(fit2_gse139929)
MODELvsCTRL <- topTable(fit2_gse139929, coef = "MODELvsCTRL", n=nrow(fit2_gse139929), adjust="fdr")
MODELvsCTRL <- data.frame(ID=rownames(MODELvsCTRL), Name=genecodev32_Annot[rownames(MODELvsCTRL),3], MODELvsCTRL)
write.table(MODELvsCTRL[MODELvsCTRL[,2] %in% c(toupper(genes_sel),"NIBAN3"),], file="GSE139929_DEGsel.txt",row.names=F, sep="\t")

cor <- apply(dataMatCounts_norm139929$E[rownames(MODELvsCTRL)[MODELvsCTRL[,2] %in% c(toupper(genes_sel),"NIBAN3")],], 1, FUN=function(x) cor(as.numeric(paste(dataMatCounts_norm139929$E["ENSG00000167483.18",])),as.numeric(paste(x))))
p <- apply(dataMatCounts_norm139929$E[rownames(MODELvsCTRL)[MODELvsCTRL[,2] %in% c(toupper(genes_sel),"NIBAN3")],], 1, FUN=function(x) cor.test(as.numeric(paste(dataMatCounts_norm139929$E["ENSG00000167483.18",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")
cor <- apply(dataMatCounts_norm139929$E[rownames(MODELvsCTRL)[MODELvsCTRL[,2] %in% c(toupper(genes_sel),"NIBAN3")],], 1, FUN=function(x) cor(as.numeric(paste(dataMatCounts_norm139929$E["ENSG00000253304.2",])),as.numeric(paste(x))))
p <- apply(dataMatCounts_norm139929$E[rownames(MODELvsCTRL)[MODELvsCTRL[,2] %in% c(toupper(genes_sel),"NIBAN3")],], 1, FUN=function(x) cor.test(as.numeric(paste(dataMatCounts_norm139929$E["ENSG00000253304.2",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")
gse139929_cor_res <- data.frame(ID=names(cor),gene=MODELvsCTRL[names(cor),2],cor=cor, p=p, fdr=fdr)
write.table(gse139929_cor_res, file="GSE139929_cor.txt", row.names=F, sep="\t", quote=F)
cor <- apply(dataMatCounts_norm139929$E[rownames(MODELvsCTRL)[MODELvsCTRL[,2] %in% c(toupper(genes_sel),"NIBAN3")],], 1, FUN=function(x) cor(as.numeric(paste(dataMatCounts_norm139929$E["ENSG00000187808.5",])),as.numeric(paste(x))))
p <- apply(dataMatCounts_norm139929$E[rownames(MODELvsCTRL)[MODELvsCTRL[,2] %in% c(toupper(genes_sel),"NIBAN3")],], 1, FUN=function(x) cor.test(as.numeric(paste(dataMatCounts_norm139929$E["ENSG00000187808.5",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")

pdf("GSE139929_TMEM200B_cor_Mar20.pdf", width = 12, height = 8)
{ 
	cor_plot <- data.frame("TMEM200B_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000253304.2",]),"TMEM200B_log2exp_2"=as.numeric(dataMatCounts_norm139929$E["ENSG00000253304.2",]))
	p1=ggplot(cor_plot, aes(x=TMEM200B_log2exp, y=TMEM200B_log2exp_2)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("TMEM200B: cor=", round(gse139929_cor_res$cor[which(gse139929_cor_res[,2]== "TMEM200B")],3)," - p=", round(gse139929_cor_res$p[which(gse139929_cor_res[,2] == "TMEM200B")],4),sep=""), x=max(cor_plot[,"TMEM200B_log2exp"]), y=max(cor_plot[,"TMEM200B_log2exp_2"]), hjust=1.8, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("TMEM200B_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000253304.2",]),"NIBAN3_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000167483.18",]))
	p2=ggplot(cor_plot, aes(x=TMEM200B_log2exp, y=NIBAN3_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("NIBAN3: cor=", round(gse139929_cor_res$cor[which(gse139929_cor_res[,2]== "NIBAN3")],3)," - p=", round(gse139929_cor_res$p[which(gse139929_cor_res[,2] == "NIBAN3")],4),sep=""), x=max(cor_plot[,"TMEM200B_log2exp"]), y=max(cor_plot[,"NIBAN3_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("TMEM200B_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000253304.2",]),"SOWAHD_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000187808.5",]))
	p3=ggplot(cor_plot, aes(x=TMEM200B_log2exp, y=SOWAHD_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("SOWAHD: cor=", round(gse139929_cor_res$cor[which(gse139929_cor_res[,2]== "SOWAHD")],3)," - p=", round(gse139929_cor_res$p[which(gse139929_cor_res[,2] == "SOWAHD")],4),sep=""), x=max(cor_plot[,"TMEM200B_log2exp"]), y=max(cor_plot[,"SOWAHD_log2exp"]), hjust=1.2, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("TMEM200B_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000253304.2",]),"WAS_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000015285.10",]))
	p4=ggplot(cor_plot, aes(x=TMEM200B_log2exp, y=WAS_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("WAS: cor=", round(gse139929_cor_res$cor[which(gse139929_cor_res[,2]== "WAS")],3)," - p=", round(gse139929_cor_res$p[which(gse139929_cor_res[,2] == "WAS")],4),sep=""), x=max(cor_plot[,"TMEM200B_log2exp"]), y=max(cor_plot[,"WAS_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("TMEM200B_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000253304.2",]),"VAV1_log2exp"=as.numeric(dataMatCounts_norm139929$E["ENSG00000141968.8",]))
	p5=ggplot(cor_plot, aes(x=TMEM200B_log2exp, y=VAV1_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("VAV1: cor=", round(gse139929_cor_res$cor[which(gse139929_cor_res[,2]== "VAV1")],3)," - p=", round(gse139929_cor_res$p[which(gse139929_cor_res[,2] == "VAV1")],4),sep=""), x=max(cor_plot[,"TMEM200B_log2exp"]), y=max(cor_plot[,"VAV1_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	multiplot(p1,p2,p3,p4,p5,cols=3)
}
dev.off()

#---------------------------------------------------
## GSE83452_NASHvsTreatNASH_n231NotPaired_Arrays_Human

gse83452 <- getGEO(filename="GSE83452/CEL/GSE83452_family.soft.gz")
charSamples83452 <- list(length(GSMList(gse83452)))
charDescription83452 <- list(length(GSMList(gse83452)))
charTitle83452 <- list(length(GSMList(gse83452)))
for (i in 1:length(GSMList(gse83452))) {                                
	charSamples83452[i] <- list(Meta(GSMList(gse83452)[[i]])$characteristics_ch1)
}
for (i in 1:length(GSMList(gse83452))) {                                
	charDescription83452[i] <- list(Meta(GSMList(gse83452)[[i]])$geo_accession)
}
for (i in 1:length(GSMList(gse83452))) {                                
	charTitle83452[i] <- list(Meta(GSMList(gse83452)[[i]])$title)
}
dataClinic83452 <- data.frame(GEO=unlist(charDescription83452),Sample=unlist(charTitle83452), Name=sapply(strsplit(sapply(charSamples83452,FUN=function(x) x[1]),": "), FUN=function(x) x[2]), LiverStatus=sapply(strsplit(sapply(charSamples83452,FUN=function(x) x[2]),": "), FUN=function(x) x[2]), Intervention=sapply(strsplit(sapply(charSamples83452,FUN=function(x) x[3]),": "), FUN=function(x) x[2]), Time=sapply(strsplit(sapply(charSamples83452,FUN=function(x) x[4]),": "), FUN=function(x) x[2]), Age=sapply(strsplit(sapply(charSamples83452,FUN=function(x) x[5]),": "), FUN=function(x) x[2]), Gender=sapply(strsplit(sapply(charSamples83452,FUN=function(x) x[6]),": "), FUN=function(x) x[2]), ScanDate=sapply(strsplit(sapply(charSamples83452,FUN=function(x) x[7]),": "), FUN=function(x) x[2]), Tissue=sapply(strsplit(sapply(charSamples83452,FUN=function(x) x[8]),": "), FUN=function(x) x[2]))
rownames(dataClinic83452) <- dataClinic83452[,1]
# NASH=126 (55_BS_baseline, 49_Diet_baseline, 1_BS_followUp, 21_Diet_followUp)
# noNASH=98 (28_BS_baseline, 16_Diet_baseline, 38_BS_followUp, 16_Diet_followUp)

gse83452_cel <- list.files(path="GSE83452/CEL",pattern=".CEL.gz", full.names=T)
gse83452_data <- read.celfiles(gse83452_cel)
gse83452_eset_rma <- rma(gse83452_data)
gse83452_mat_rma <- exprs(gse83452_eset_rma)
gene83452 <- lookUp(rownames(gse83452_mat_rma), "hugene20", "SYMBOL");

tipo83452_1 <- factor(ifelse(substr(colnames(gse83452_mat_rma), 1, 10) %in% rownames(dataClinic83452)[dataClinic83452$LiverStatus=="NASH"], "NASH","noNASH"))
tipo83452_2 <- factor(ifelse(substr(colnames(gse83452_mat_rma), 1, 10) %in% rownames(dataClinic83452)[dataClinic83452$Intervention=="BS"], "BS","Diet"))
tipo83452_3 <- factor(ifelse(substr(colnames(gse83452_mat_rma), 1, 10) %in% rownames(dataClinic83452)[dataClinic83452$Time=="baseline"], "baseline","followUp"))
tipo83452 <- paste(tipo83452_1,tipo83452_2,tipo83452_3,sep="_")
gse83452_mat_rma_filter <- gse83452_mat_rma[apply(gse83452_mat_rma[,tipo83452=="NASH_BS_baseline"],1,FUN=function(x) sum(x>5))>=28 | apply(gse83452_mat_rma[,tipo83452=="NASH_Diet_baseline"],1,FUN=function(x) sum(x>5))>=25 | apply(gse83452_mat_rma[,tipo83452=="NASH_Diet_followup"],1,FUN=function(x) sum(x>5))>=11 | apply(gse83452_mat_rma[,tipo83452=="noNASH_BS_baseline"],1,FUN=function(x) sum(x>5))>=14 | apply(gse83452_mat_rma[,tipo83452=="noNASH_BS_followup"],1,FUN=function(x) sum(x>5))>=19 | apply(gse83452_mat_rma[,tipo83452=="noNASH_Diet_baseline"],1,FUN=function(x) sum(x>5))>=8 | apply(gse83452_mat_rma[,tipo83452=="noNASH_Diet_followup"],1,FUN=function(x) sum(x>5))>=8,]

load("~/dato-activo/00_References/AffyAnnot/HuGene20stv1_na36_hg19.rda")
gse83452_mat_rma_filter_annot <- merge(ann,data.frame(Probeset=rownames(gse83452_mat_rma_filter), gse83452_mat_rma_filter), by.x = 1, by.y = 1)

# QC
acol = sample(brewer.pal(8, "Dark2"), ncol(gse83452_mat_rma), replace = (8 < ncol(gse83452_mat_rma)))
xlim_r = c(min(na.omit(exprs(gse83452_data))),max(na.omit(exprs(gse83452_data)))) 
xlim = c(min(na.omit(gse83452_mat_rma)),max(na.omit(gse83452_mat_rma)))
clust.euclid.average <- hclust(dist(t(exprs(gse83452_data))),method="average")
clust.euclid.average2 <- hclust(dist(t(gse83452_mat_rma)),method="average")
pdf(file="GSE83452_QC_Array_Mar20.pdf", width = 20, height = 20)
boxplot(exprs(gse83452_data),names=colnames(gse83452_eset_rma), which="all", transfo=log2, cex.axis=0.7, las=2)
multi("density", exprs(gse83452_data), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(gse83452_mat_rma, cex.axis=0.7, las=2)
multi("density", gse83452_mat_rma, xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)
dev.off()

a <- detOutliers(gse83452_mat_rma, tipo83452, "GSE83452_Array_var_Mar20.pdf", 20, 20)
a <- detOutliers(gse83452_mat_rma, tipo83452_1, "GSE83452_Array_var2_Mar20.pdf", 20, 20)
a <- detOutliers(gse83452_mat_rma, tipo83452_3, "GSE83452_Array_var2b_Mar20.pdf", 20, 20)

x="PC1"
y="PC2"
fit <- prcomp(t(na.omit(as.matrix(dataMatCounts_norm139929$E))), scale=T)
pcData <- data.frame(fit$x)
ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=tipo139929)) + geom_text(aes(label=colnames(dataMatCounts_norm139929$E)), hjust=0, vjust=0, nudge_x=-8, nudge_y=2) + theme_bw() + theme(legend.title=element_blank())
pdf(file = "GSE139929_PCA_Mar20.pdf", width = 10, height = 10, colormodel = "rgb")
ggp
dev.off()

design_gse83452 <- model.matrix(~0+tipo83452)
rownames(design_gse83452) <- colnames(gse83452_mat_rma_filter)
colnames(design_gse83452) <- levels(tipo83452)
cont.matrix_gse83452 <- makeContrasts(NASHdietFUvsNASHdiet = NASH_Diet_followUp - NASH_Diet_baseline,
	noNASHdietFUvsNoNASHdiet = noNASH_Diet_followUp - noNASH_Diet_baseline,
	noNASHbsFUvsNoNASHbs = noNASH_BS_followUp - noNASH_BS_baseline,
	NASHdietFUvsNoNASHdietFU = NASH_Diet_followUp - noNASH_Diet_followUp,
	NASHdietvsNoNASHdiet = NASH_Diet_baseline - noNASH_Diet_baseline,
	NASHbsvsNoNASHbs = NASH_BS_baseline - noNASH_BS_baseline,
	levels=design_gse83452)
fit_gse83452 <- lmFit(gse83452_mat_rma_filter, design_gse83452)
fit2_gse83452 <- contrasts.fit(fit_gse83452, cont.matrix_gse83452)
fit2_gse83452 <- eBayes(fit2_gse83452)
NASHdietFUvsNASHdiet <- topTable(fit2_gse83452, coef = "NASHdietFUvsNASHdiet", n=nrow(fit2_gse83452), adjust="fdr")
NASHdietFUvsNASHdiet <- merge(ann,data.frame(Probeset=rownames(NASHdietFUvsNASHdiet), NASHdietFUvsNASHdiet), by.x = 1, by.y = 1)
NASHdietFUvsNASHdiet <- NASHdietFUvsNASHdiet[order(NASHdietFUvsNASHdiet$P.Value),]
NASHdietFUvsNASHdiet[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel),]
rownames(NASHdietFUvsNASHdiet) <- NASHdietFUvsNASHdiet[,1]
noNASHdietFUvsNoNASHdiet <- topTable(fit2_gse83452, coef = "noNASHdietFUvsNoNASHdiet", n=nrow(fit2_gse83452), adjust="fdr")
noNASHdietFUvsNoNASHdiet <- merge(ann,data.frame(Probeset=rownames(noNASHdietFUvsNoNASHdiet), noNASHdietFUvsNoNASHdiet), by.x = 1, by.y = 1)
noNASHdietFUvsNoNASHdiet <- noNASHdietFUvsNoNASHdiet[order(noNASHdietFUvsNoNASHdiet$P.Value),]
noNASHdietFUvsNoNASHdiet[noNASHdietFUvsNoNASHdiet[,3] %in% toupper(genes_sel),]
rownames(noNASHdietFUvsNoNASHdiet) <- noNASHdietFUvsNoNASHdiet[,1]
noNASHbsFUvsNoNASHbs <- topTable(fit2_gse83452, coef = "noNASHbsFUvsNoNASHbs", n=nrow(fit2_gse83452), adjust="fdr")
noNASHbsFUvsNoNASHbs <- merge(ann,data.frame(Probeset=rownames(noNASHbsFUvsNoNASHbs), noNASHbsFUvsNoNASHbs), by.x = 1, by.y = 1)
noNASHbsFUvsNoNASHbs <- noNASHbsFUvsNoNASHbs[order(noNASHbsFUvsNoNASHbs$P.Value),]
noNASHbsFUvsNoNASHbs[noNASHbsFUvsNoNASHbs[,3] %in% toupper(genes_sel),]
rownames(noNASHbsFUvsNoNASHbs) <- noNASHbsFUvsNoNASHbs[,1]
NASHdietFUvsNoNASHdietFU <- topTable(fit2_gse83452, coef = "NASHdietFUvsNoNASHdietFU", n=nrow(fit2_gse83452), adjust="fdr")
NASHdietFUvsNoNASHdietFU <- merge(ann,data.frame(Probeset=rownames(NASHdietFUvsNoNASHdietFU), NASHdietFUvsNoNASHdietFU), by.x = 1, by.y = 1)
NASHdietFUvsNoNASHdietFU <- NASHdietFUvsNoNASHdietFU[order(NASHdietFUvsNoNASHdietFU$P.Value),]
NASHdietFUvsNoNASHdietFU[NASHdietFUvsNoNASHdietFU[,3] %in% toupper(genes_sel),]
rownames(NASHdietFUvsNoNASHdietFU) <- NASHdietFUvsNoNASHdietFU[,1]
NASHdietvsNoNASHdiet <- topTable(fit2_gse83452, coef = "NASHdietvsNoNASHdiet", n=nrow(fit2_gse83452), adjust="fdr")
NASHdietvsNoNASHdiet <- merge(ann,data.frame(Probeset=rownames(NASHdietvsNoNASHdiet), NASHdietvsNoNASHdiet), by.x = 1, by.y = 1)
NASHdietvsNoNASHdiet <- NASHdietvsNoNASHdiet[order(NASHdietvsNoNASHdiet$P.Value),]
NASHdietvsNoNASHdiet[NASHdietvsNoNASHdiet[,3] %in% toupper(genes_sel),]
rownames(NASHdietvsNoNASHdiet) <- NASHdietvsNoNASHdiet[,1]
NASHbsvsNoNASHbs <- topTable(fit2_gse83452, coef = "NASHbsvsNoNASHbs", n=nrow(fit2_gse83452), adjust="fdr")
NASHbsvsNoNASHbs <- merge(ann,data.frame(Probeset=rownames(NASHbsvsNoNASHbs), NASHbsvsNoNASHbs), by.x = 1, by.y = 1)
NASHbsvsNoNASHbs <- NASHbsvsNoNASHbs[order(NASHbsvsNoNASHbs$P.Value),]
NASHbsvsNoNASHbs[NASHbsvsNoNASHbs[,3] %in% toupper(genes_sel),]
rownames(NASHbsvsNoNASHbs) <- NASHbsvsNoNASHbs[,1]

design_gse83452_1 <- model.matrix(~0+tipo83452_1)
rownames(design_gse83452_1) <- colnames(gse83452_mat_rma_filter)
colnames(design_gse83452_1) <- levels(tipo83452_1)
cont.matrix_gse83452_1 <- makeContrasts(NASHvsNoNASH = NASH - noNASH,
	levels=design_gse83452_1)
fit_gse83452_1 <- lmFit(gse83452_mat_rma_filter, design_gse83452_1)
fit2_gse83452_1 <- contrasts.fit(fit_gse83452_1, cont.matrix_gse83452_1)
fit2_gse83452_1 <- eBayes(fit2_gse83452_1)
NASHvsNoNASH <- topTable(fit2_gse83452_1, coef = "NASHvsNoNASH", n=nrow(fit2_gse83452_1), adjust="fdr")
NASHvsNoNASH <- merge(ann,data.frame(Probeset=rownames(NASHvsNoNASH), NASHvsNoNASH), by.x = 1, by.y = 1)
NASHvsNoNASH <- NASHvsNoNASH[order(NASHvsNoNASH$P.Value),]
NASHvsNoNASH[NASHvsNoNASH[,3] %in% toupper(genes_sel),]
rownames(NASHvsNoNASH) <- NASHvsNoNASH[,1]

design_gse83452_3 <- model.matrix(~0+tipo83452_3)
rownames(design_gse83452_3) <- colnames(gse83452_mat_rma_filter)
colnames(design_gse83452_3) <- levels(tipo83452_3)
cont.matrix_gse83452_3 <- makeContrasts(FollowUPvsBaseline = followUp - baseline,
	levels=design_gse83452_3)
fit_gse83452_3 <- lmFit(gse83452_mat_rma_filter, design_gse83452_3)
fit2_gse83452_3 <- contrasts.fit(fit_gse83452_3, cont.matrix_gse83452_3)
fit2_gse83452_3 <- eBayes(fit2_gse83452_3)
FollowUPvsBaseline <- topTable(fit2_gse83452_3, coef = "FollowUPvsBaseline", n=nrow(fit2_gse83452_3), adjust="fdr")
FollowUPvsBaseline <- merge(ann,data.frame(Probeset=rownames(FollowUPvsBaseline), FollowUPvsBaseline), by.x = 1, by.y = 1)
FollowUPvsBaseline <- FollowUPvsBaseline[order(FollowUPvsBaseline$P.Value),]
FollowUPvsBaseline[FollowUPvsBaseline[,3] %in% c(paste(metabolism_uPE_info[,3]),"FAM129C"),]
rownames(FollowUPvsBaseline) <- FollowUPvsBaseline[,1]
write.table(FollowUPvsBaseline, "FollowUPvsBaseline_GSE83452.txt", row.names=FALSE, sep="\t", quote=F)

GSE83452_bmat <- cbind(NASHdietFUvsNASHdiet[rownames(NASHdietFUvsNASHdiet)[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel)],c(1,3,8,11,12)], noNASHdietFUvsNoNASHdiet[rownames(NASHdietFUvsNASHdiet)[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel)],c(8,11,12)],noNASHbsFUvsNoNASHbs[rownames(NASHdietFUvsNASHdiet)[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel)],c(8,11,12)], NASHdietFUvsNoNASHdietFU[rownames(NASHdietFUvsNASHdiet)[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel)],c(8,11,12)], NASHdietvsNoNASHdiet[rownames(NASHdietFUvsNASHdiet)[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel)],c(8,11,12)], NASHbsvsNoNASHbs[rownames(NASHdietFUvsNASHdiet)[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel)],c(8,11,12)], NASHvsNoNASH[rownames(NASHdietFUvsNASHdiet)[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel)],c(8,11,12)], FollowUPvsBaseline[rownames(NASHdietFUvsNASHdiet)[NASHdietFUvsNASHdiet[,3] %in% toupper(genes_sel)],c(8,11,12)])
colnames(GSE83452_bmat)[-c(1:2)] <- paste(c(rep("NASHdietFUvsNASHdiet",3), rep("noNASHdietFUvsNoNASHdiet",3), rep("noNASHbsFUvsNoNASHbs",3),rep("NASHdietFUvsNoNASHdietFU",3),rep("NASHdietvsNoNASHdiet",3), rep("NASHbsvsNoNASHbs",3), rep("NASHvsNoNASH",3), rep("FollowUPvsBaseline",3)),colnames(GSE83452_bmat)[-c(1:2)], sep=".")
write.table(GSE83452_bmat, file="GSE83452_DEGsel_allSamples.txt",row.names=F, sep="\t")

cor <- apply(gse83452_mat_rma_filter[paste(FollowUPvsBaseline[FollowUPvsBaseline[,3] %in% c(toupper(genes_sel),"NIBAN3"),1]),], 1, FUN=function(x) cor(as.numeric(paste(gse83452_mat_rma_filter["16672082",])),as.numeric(paste(x))))
p <- apply(gse83452_mat_rma_filter[paste(FollowUPvsBaseline[FollowUPvsBaseline[,3] %in% c(toupper(genes_sel),"NIBAN3"),1]),], 1, FUN=function(x) cor.test(as.numeric(paste(gse83452_mat_rma_filter["16672082",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")
gse83452followUp_cor_res <- data.frame(ID=names(cor),gene=NASHvsNoNASH[names(cor),3],cor=cor, p=p, fdr=fdr)
write.table(gse83452followUp_cor_res, file="GSE83452_TTC24_allSamples_cor.txt", row.names=F, sep="\t", quote=F)
cor <- apply(gse83452_mat_rma_filter[paste(FollowUPvsBaseline[FollowUPvsBaseline[,3] %in% c(toupper(genes_sel),"NIBAN3"),1]),], 1, FUN=function(x) cor(as.numeric(paste(gse83452_mat_rma_filter["16859562",])),as.numeric(paste(x))))
p <- apply(gse83452_mat_rma_filter[paste(FollowUPvsBaseline[FollowUPvsBaseline[,3] %in% c(toupper(genes_sel),"NIBAN3"),1]),], 1, FUN=function(x) cor.test(as.numeric(paste(gse83452_mat_rma_filter["16859562",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")
gse83452followUp_cor_res2 <- data.frame(ID=names(cor),gene=NASHvsNoNASH[names(cor),3],cor=cor, p=p, fdr=fdr)
write.table(gse83452followUp_cor_res2, file="GSE83452_FAM129C_allSamples_cor.txt", row.names=F, sep="\t", quote=F)
pdf("GSE83452_FAM129C_cor_Allsamples_Mar20.pdf", width = 12, height = 8)
{ 
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(gse83452_mat_rma_filter["16859562",]),"FAM129C_log2exp_2"=as.numeric(gse83452_mat_rma_filter["16859562",]))
	p1=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=FAM129C_log2exp_2)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("FAM129C: cor=", round(gse83452followUp_cor_res2$cor[which(gse83452followUp_cor_res2[,2]== "FAM129C")],3)," - p=", round(gse83452followUp_cor_res2$p[which(gse83452followUp_cor_res2[,2] == "FAM129C")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"FAM129C_log2exp_2"]), hjust=1.8, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(gse83452_mat_rma_filter["16859562",]),"TTC24_log2exp"=as.numeric(gse83452_mat_rma_filter["16672082",]))
	p2=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=TTC24_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("TTC24: cor=", round(gse83452followUp_cor_res2$cor[which(gse83452followUp_cor_res2[,2]== "TTC24")],3)," - p=", round(gse83452followUp_cor_res2$p[which(gse83452followUp_cor_res2[,2] == "TTC24")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"TTC24_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(gse83452_mat_rma_filter["16859562",]),"PLEK_log2exp"=as.numeric(gse83452_mat_rma_filter["16880942",]))
	p3=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=PLEK_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("PLEK: cor=", round(gse83452followUp_cor_res2$cor[which(gse83452followUp_cor_res2[,2]== "PLEK")],3)," - p=", round(gse83452followUp_cor_res2$p[which(gse83452followUp_cor_res2[,2] == "PLEK")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"PLEK_log2exp"]), hjust=1.2, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(gse83452_mat_rma_filter["16859562",]),"INPP5D_log2exp"=as.numeric(gse83452_mat_rma_filter["16892446",]))
	p4=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=INPP5D_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("INPP5D: cor=", round(gse83452followUp_cor_res2$cor[which(gse83452followUp_cor_res2[,2]== "INPP5D")],3)," - p=", round(gse83452followUp_cor_res2$p[which(gse83452followUp_cor_res2[,2] == "INPP5D")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"INPP5D_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	cor_plot <- data.frame("FAM129C_log2exp"=as.numeric(gse83452_mat_rma_filter["16859562",]),"VAV1_log2exp"=as.numeric(gse83452_mat_rma_filter["16857490",]))
	p5=ggplot(cor_plot, aes(x=FAM129C_log2exp, y=VAV1_log2exp)) + geom_point() + geom_smooth(method=lm, se=FALSE, color="black") + theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) + annotate("text", label=paste("VAV1: cor=", round(gse83452followUp_cor_res2$cor[which(gse83452followUp_cor_res2[,2]== "VAV1")],3)," - p=", round(gse83452followUp_cor_res2$p[which(gse83452followUp_cor_res2[,2] == "VAV1")],4),sep=""), x=max(cor_plot[,"FAM129C_log2exp"]), y=max(cor_plot[,"VAV1_log2exp"]), hjust=1.6, size=4, vjust=0) + theme_bw()
	multiplot(p1,p2,p3,p4,p5,cols=3)
}
dev.off()

# busco gene-set de genes correlados para FAM129
NASHdietFUvsNASHdiet[NASHdietFUvsNASHdiet[,3] %in% paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]) & NASHdietFUvsNASHdiet$P.Value<0.05,]
GSE83452_FAM129C_phyper <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(NASHdietFUvsNASHdiet[,3])), paste(NASHdietFUvsNASHdiet[NASHdietFUvsNASHdiet$P.Value<0.05,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(NASHdietFUvsNASHdiet[,3]))), length(unique(NASHdietFUvsNASHdiet[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(NASHdietFUvsNASHdiet[,3]))), length(unique(NASHdietFUvsNASHdiet[NASHdietFUvsNASHdiet$P.Value<0.05,3])), lower.tail=F) # 0.17
GSE83452_FAM129C_phyper2 <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(NASHdietFUvsNASHdiet[,3])), paste(NASHdietFUvsNASHdiet[NASHdietFUvsNASHdiet$P.Value<0.01,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(NASHdietFUvsNASHdiet[,3]))), length(unique(NASHdietFUvsNASHdiet[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(NASHdietFUvsNASHdiet[,3]))), length(unique(NASHdietFUvsNASHdiet[NASHdietFUvsNASHdiet$P.Value<0.01,3])), lower.tail=F) # 0.68
GSE83452_FAM129C_phyperB <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3])), paste(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.05,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.05,3])), lower.tail=F) # 0.17
GSE83452_FAM129C_phyper2B <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3])), paste(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.01,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.01,3])), lower.tail=F) # 0.68
GSE83452_FAM129C_phyper3B <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3])), paste(FollowUPvsBaseline[FollowUPvsBaseline$adj.P.Val<0.05,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[FollowUPvsBaseline$adj.P.Val<0.05,3])), lower.tail=F) # 0.68


pdf("GSE83452_FAM129Cenrichment.pdf",16,8)
par(mfrow=c(1,2))
compare2List(paste(unique(NASHdietFUvsNASHdiet[NASHdietFUvsNASHdiet$P.Value<0.05,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(NASHdietFUvsNASHdiet[,3])), "NASHdietFUvsNASHdiet_p05","FAM129Ccor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE83452_FAM129C_phyper,3),sep=""))
compare2List(paste(unique(NASHdietFUvsNASHdiet[NASHdietFUvsNASHdiet$P.Value<0.01,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(NASHdietFUvsNASHdiet[,3])), "NASHdietFUvsNASHdiet_p01","FAM129Ccor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE83452_FAM129C_phyper2,3),sep=""))
compare2List(paste(unique(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.05,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3])), "FollowUPvsBaseline_p05","FAM129Ccor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE83452_FAM129C_phyperB,5),sep=""))
compare2List(paste(unique(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.01,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(FollowUPvsBaseline[,3])), "FollowUPvsBaseline_p01","FAM129Ccor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE83452_FAM129C_phyper2B,7),sep=""))
dev.off()


# busco gene-set de genes correlados para TTC24
TTC24_corTCGA <- rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,"ENSG00000187862"]>0]
TTC24_corCCLE <- rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,"ENSG00000187862"]>0]
TTC24_corGTEX <- rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,"ENSG00000187862"]>0]
TTC24_cor <- intersect(TTC24_corTCGA,intersect(TTC24_corCCLE,TTC24_corGTEX))

genecodev28_all <- read.table("~/dato-activo/00_References/gencode_gtf/gencode.v28.annotation.gtf", skip = 5, header = FALSE, sep = "\t")
genecodev28 <- genecodev28_all[genecodev28_all$V3 == "gene", ]
colnames(genecodev28) <- c("chr","DB","Type", "start","end","","strand","","Description")
genecodev28_tmp <- apply(as.data.frame(genecodev28[,9]), 1, parseENCODE_eli)   
genecodev28_tmp <- t(genecodev28_tmp)
colnames(genecodev28_tmp) <- c("gene_id", "gene_type", "gene_name", "level", "havana_gene")
genecodev28_Annot <- data.frame(genecodev28_tmp, genecodev28[,-c(6,8:9)])
rownames(genecodev28_Annot) <- genecodev28_Annot[,1]
genecodev28_Annot$ENSG <- substr(genecodev28_Annot[,1],1,15)


GSE83452_TTC24_phyper <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(noNASHbsFUvsNoNASHbs[,3])), paste(noNASHbsFUvsNoNASHbs[noNASHbsFUvsNoNASHbs$P.Value<0.05,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(noNASHbsFUvsNoNASHbs[,3]))), length(unique(noNASHbsFUvsNoNASHbs[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(noNASHbsFUvsNoNASHbs[,3]))), length(unique(noNASHbsFUvsNoNASHbs[noNASHbsFUvsNoNASHbs$P.Value<0.05,3])), lower.tail=F) # 0.00035
GSE83452_TTC24_phyper2 <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(noNASHbsFUvsNoNASHbs[,3])), paste(noNASHbsFUvsNoNASHbs[noNASHbsFUvsNoNASHbs$P.Value<0.01,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(noNASHbsFUvsNoNASHbs[,3]))), length(unique(noNASHbsFUvsNoNASHbs[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(noNASHbsFUvsNoNASHbs[,3]))), length(unique(noNASHbsFUvsNoNASHbs[noNASHbsFUvsNoNASHbs$P.Value<0.01,3])), lower.tail=F) # 1.14E-05

GSE83452_TTC24_phyperB <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3])), paste(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.05,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.05,3])), lower.tail=F) # 0.00035 PAPER: 0.004702952
GSE83452_TTC24_phyper2B <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3])), paste(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.01,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.01,3])), lower.tail=F) # 1.14E-05 PAPER: 0.002260
GSE83452_TTC24_phyper3B <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3])), paste(FollowUPvsBaseline[FollowUPvsBaseline$adj.P.Val<0.05,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[FollowUPvsBaseline$adj.P.Val<0.05,3])), lower.tail=F) # 1.14E-05 PAPER: 0.002090


pdf("GSE83452_TTC24enrichment.pdf",16,8)
par(mfrow=c(1,2))
compare2List(paste(unique(noNASHbsFUvsNoNASHbs[noNASHbsFUvsNoNASHbs$P.Value<0.05,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% FAM129C_cor,3]),paste(noNASHbsFUvsNoNASHbs[,3])), "noNASHbsFUvsNoNASHbs_p05","TTC24cor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE83452_TTC24_phyper,5),sep=""))
compare2List(paste(unique(noNASHbsFUvsNoNASHbs[noNASHbsFUvsNoNASHbs$P.Value<0.01,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(noNASHbsFUvsNoNASHbs[,3])), "noNASHbsFUvsNoNASHbs_p01","TTC24cor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE83452_TTC24_phyper2,7),sep=""))
compare2List(paste(unique(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.05,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3])), "FollowUPvsBaseline_p05","TTC24cor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE83452_TTC24_phyperB,5),sep=""))
compare2List(paste(unique(FollowUPvsBaseline[FollowUPvsBaseline$P.Value<0.01,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3])), "FollowUPvsBaseline_p01","TTC24cor_GTEX.TCGA.CCLE",paste("phyper=",round(GSE83452_TTC24_phyper2B,7),sep=""))
dev.off()

# paper Figure
pdf("GSE83452_EnrichmentStudy.pdf",4,4)
compare2List(paste(unique(FollowUPvsBaseline[FollowUPvsBaseline$adj.P.Val<0.05,3])), intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% TTC24_cor,3]),paste(FollowUPvsBaseline[,3])), "Obesity treatment","TTC24 correlated genes",paste("p=",round(GSE83452_TTC24_phyper3B,5),sep=""))
dev.off()

# busco gene-set de genes correlados para CPED1
CPED1_corTCGA <- rownames(TCGA_uPExPE1)[TCGA_uPExPE1[,"ENSG00000106034"]>0]
CPED1_corCCLE <- rownames(CCLE_uPExPE1)[CCLE_uPExPE1[,"ENSG00000106034"]>0]
CPED1_corGTEX <- rownames(GTEX_uPExPE1)[GTEX_uPExPE1[,"ENSG00000106034"]>0]
CPED1_cor <- intersect(CPED1_corTCGA,intersect(CPED1_corCCLE,CPED1_corGTEX))

GSE83452_CPED1_phyper3B <- phyper(length(intersect(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% CPED1_cor,3]),paste(FollowUPvsBaseline[,3])), paste(FollowUPvsBaseline[FollowUPvsBaseline$adj.P.Val<0.05,3]))), length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% CPED1_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[,3]))-length(intersect(paste(genecodev28_Annot[genecodev28_Annot$ENSG %in% CPED1_cor,3]),paste(FollowUPvsBaseline[,3]))), length(unique(FollowUPvsBaseline[FollowUPvsBaseline$adj.P.Val<0.05,3])), lower.tail=F) # 1.14E-05 PAPER: 0.7991


## Only NASH diet samples
gse83452_data2 <- read.celfiles(gse83452_cel[tipo83452 %in% c("NASH_Diet_baseline","NASH_Diet_followUp")])
gse83452_eset_rma2 <- rma(gse83452_data2)
gse83452_mat_rma2 <- exprs(gse83452_eset_rma2)
tipo83452_3b <- factor(ifelse(substr(colnames(gse83452_mat_rma2), 1, 10) %in% rownames(dataClinic83452)[dataClinic83452$Time=="baseline"], "baseline","followUp"))
gse83452_mat_rma_filter2 <- gse83452_mat_rma2[apply(gse83452_mat_rma[,tipo83452=="NASH_Diet_baseline"],1,FUN=function(x) sum(x>5))>=25 | apply(gse83452_mat_rma[,tipo83452=="NASH_Diet_followup"],1,FUN=function(x) sum(x>5))>=11,]

design_gse83452b <- model.matrix(~0+tipo83452_3b)
rownames(design_gse83452b) <- colnames(gse83452_mat_rma_filter2)
colnames(design_gse83452b) <- levels(tipo83452_3b)
cont.matrix_gse83452b <- makeContrasts(NASHdietFUvsNASHdiet_b = followUp - baseline,
		levels=design_gse83452b)
fit_gse83452b <- lmFit(gse83452_mat_rma_filter2, design_gse83452b)
fit2_gse83452b <- contrasts.fit(fit_gse83452b, cont.matrix_gse83452b)
fit2_gse83452b <- eBayes(fit2_gse83452b)
NASHdietFUvsNASHdiet_b <- topTable(fit2_gse83452b, coef = "NASHdietFUvsNASHdiet_b", n=nrow(fit2_gse83452b), adjust="fdr")
NASHdietFUvsNASHdiet_b <- merge(ann,data.frame(Probeset=rownames(NASHdietFUvsNASHdiet_b), NASHdietFUvsNASHdiet_b), by.x = 1, by.y = 1)
NASHdietFUvsNASHdiet_b <- NASHdietFUvsNASHdiet_b[order(NASHdietFUvsNASHdiet_b$P.Value),]
NASHdietFUvsNASHdiet_b[NASHdietFUvsNASHdiet_b[,3] %in% toupper(genes_sel),]
write.table(NASHdietFUvsNASHdiet_b[NASHdietFUvsNASHdiet_b[,2] %in% toupper(genes_sel),], file="GSE83452_NASHdietFUvsNASHdiet_DEGsel_b.txt",row.names=F, sep="\t")
cor <- apply(gse83452_mat_rma_filter2[paste(FollowUPvsBaseline[FollowUPvsBaseline[,3] %in% c(toupper(genes_sel),"NIBAN3"),1]),], 1, FUN=function(x) cor(as.numeric(paste(gse83452_mat_rma_filter2["16859562",])),as.numeric(paste(x))))
p <- apply(gse83452_mat_rma_filter2[paste(FollowUPvsBaseline[FollowUPvsBaseline[,3] %in% c(toupper(genes_sel),"NIBAN3"),1]),], 1, FUN=function(x) cor.test(as.numeric(paste(gse83452_mat_rma_filter2["16859562",])),as.numeric(paste(x)))$p.value)
fdr <- p.adjust(p, method="fdr")
gse83452followUp_cor_res2b <- data.frame(ID=names(cor),gene=NASHvsNoNASH[names(cor),3],cor=cor, p=p, fdr=fdr)

#####################################################################################################################################################

## Nostromo (Function distance)
setwd("/home/nostromo/data/03_Analysis/eguruce/20_UPs_Abr19")

# load("EvaluationOfPR_sep20.RData")
# save.image("EvaluationOfPR_sep20.RData")

##########################
# New known PE1 
##########################

##########################
# New known PE1 x nextProt
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/CCLE/",full.names=T)
newPE_PRvsNx_CCLEBP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_CCLEBP <- tmp

  } else {

    newPE_PRvsNx_CCLEBP <- rbind(newPE_PRvsNx_CCLEBP, tmp)

  }

}
colnames(newPE_PRvsNx_CCLEBP) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_CCLEBP$Dist <- as.numeric(paste(newPE_PRvsNx_CCLEBP$Dist))
newPE_PRvsNx_CCLEBP_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_CCLEBP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/CCLE/",full.names=T)
newPE_PRvsNx_CCLEMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_CCLEMF <- tmp

  } else {

    newPE_PRvsNx_CCLEMF <- rbind(newPE_PRvsNx_CCLEMF, tmp)

  }

}
colnames(newPE_PRvsNx_CCLEMF) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_CCLEMF$Dist <- as.numeric(paste(newPE_PRvsNx_CCLEMF$Dist))
newPE_PRvsNx_CCLEMF_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_CCLEMF[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/CCLE/",full.names=T)
newPE_PRvsNx_CCLECC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_CCLECC <- tmp

  } else {

    newPE_PRvsNx_CCLECC <- rbind(newPE_PRvsNx_CCLECC, tmp)

  }

}
colnames(newPE_PRvsNx_CCLECC) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_CCLECC$Dist <- as.numeric(paste(newPE_PRvsNx_CCLECC$Dist))
newPE_PRvsNx_CCLECC_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_CCLECC[,c(1,5)], FUN=max, na.rm=T)

newPE_PRvsNx_CCLE <- rbind(newPE_PRvsNx_CCLEBP_max,newPE_PRvsNx_CCLEMF_max,newPE_PRvsNx_CCLECC_max)
newPE_PRvsNx_CCLE[,1] <- sapply(paste(newPE_PRvsNx_CCLE[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
newPE_PRvsNx_CCLE_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_CCLE, FUN=max, na.rm=T)
colnames(newPE_PRvsNx_CCLE_max)[2] <- "Dist.max"

# TCGA
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/TCGA/",full.names=T)
newPE_PRvsNx_TCGABP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_TCGABP <- tmp

  } else {

    newPE_PRvsNx_TCGABP <- rbind(newPE_PRvsNx_TCGABP, tmp)

  }

}
colnames(newPE_PRvsNx_TCGABP) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_TCGABP$Dist <- as.numeric(paste(newPE_PRvsNx_TCGABP$Dist))
newPE_PRvsNx_TCGABP_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_TCGABP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/TCGA/",full.names=T)
newPE_PRvsNx_TCGAMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_TCGAMF <- tmp

  } else {

    newPE_PRvsNx_TCGAMF <- rbind(newPE_PRvsNx_TCGAMF, tmp)

  }

}
colnames(newPE_PRvsNx_TCGAMF) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_TCGAMF$Dist <- as.numeric(paste(newPE_PRvsNx_TCGAMF$Dist))
newPE_PRvsNx_TCGAMF_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_TCGAMF[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/TCGA/",full.names=T)
newPE_PRvsNx_TCGACC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_TCGACC <- tmp

  } else {

    newPE_PRvsNx_TCGACC <- rbind(newPE_PRvsNx_TCGACC, tmp)

  }

}
colnames(newPE_PRvsNx_TCGACC) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_TCGACC$Dist <- as.numeric(paste(newPE_PRvsNx_TCGACC$Dist))
newPE_PRvsNx_TCGACC_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_TCGACC[,c(1,5)], FUN=max, na.rm=T)

newPE_PRvsNx_TCGA <- rbind(newPE_PRvsNx_TCGABP_max,newPE_PRvsNx_TCGAMF_max,newPE_PRvsNx_TCGACC_max)
newPE_PRvsNx_TCGA[,1] <- sapply(paste(newPE_PRvsNx_TCGA[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
newPE_PRvsNx_TCGA_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_TCGA, FUN=max, na.rm=T)
colnames(newPE_PRvsNx_TCGA_max)[2] <- "Dist.max"

# GTEX
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/GTEX/",full.names=T)
newPE_PRvsNx_GTEXBP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_GTEXBP <- tmp

  } else {

    newPE_PRvsNx_GTEXBP <- rbind(newPE_PRvsNx_GTEXBP, tmp)

  }

}
colnames(newPE_PRvsNx_GTEXBP) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_GTEXBP$Dist <- as.numeric(paste(newPE_PRvsNx_GTEXBP$Dist))
newPE_PRvsNx_GTEXBP_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_GTEXBP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/GTEX/",full.names=T)
newPE_PRvsNx_GTEXMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_GTEXMF <- tmp

  } else {

    newPE_PRvsNx_GTEXMF <- rbind(newPE_PRvsNx_GTEXMF, tmp)

  }

}
colnames(newPE_PRvsNx_GTEXMF) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_GTEXMF$Dist <- as.numeric(paste(newPE_PRvsNx_GTEXMF$Dist))
newPE_PRvsNx_GTEXMF_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_GTEXMF[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_NEXTPROT/GTEX/",full.names=T)
newPE_PRvsNx_GTEXCC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsNx_GTEXCC <- tmp

  } else {

    newPE_PRvsNx_GTEXCC <- rbind(newPE_PRvsNx_GTEXCC, tmp)

  }

}
colnames(newPE_PRvsNx_GTEXCC) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsNx_GTEXCC$Dist <- as.numeric(paste(newPE_PRvsNx_GTEXCC$Dist))
newPE_PRvsNx_GTEXCC_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_GTEXCC[,c(1,5)], FUN=max, na.rm=T)

newPE_PRvsNx_GTEX <- rbind(newPE_PRvsNx_GTEXBP_max,newPE_PRvsNx_GTEXMF_max,newPE_PRvsNx_GTEXCC_max)
newPE_PRvsNx_GTEX[,1] <- sapply(paste(newPE_PRvsNx_GTEX[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
newPE_PRvsNx_GTEX_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx_GTEX, FUN=max, na.rm=T)
colnames(newPE_PRvsNx_GTEX_max)[2] <- "Dist.max"

## ALL DBs and Datasets
newPE_PRvsNx <- rbind(newPE_PRvsNx_CCLE_max,newPE_PRvsNx_TCGA_max,newPE_PRvsNx_GTEX_max)
newPE_PRvsNx_max <- summaryBy(. ~ Prot, data=newPE_PRvsNx, FUN=max, na.rm=T)
colnames(newPE_PRvsNx_max)[2] <- "Dist.max"

newPE_PRvsNx_ggplot2 <- rbind(data.frame(newPE_PRvsNx_CCLEBP_max,DB="GO_BP",Experiment="CCLE"),
  data.frame(newPE_PRvsNx_CCLEMF_max,DB="GO_MF",Experiment="CCLE"),
 data.frame(newPE_PRvsNx_CCLECC_max,DB="GO_CC",Experiment="CCLE"),
 data.frame(newPE_PRvsNx_CCLE_max,DB="All",Experiment="CCLE"),
 data.frame(newPE_PRvsNx_GTEXBP_max,DB="GO_BP",Experiment="GTEX"),
 data.frame(newPE_PRvsNx_GTEXMF_max,DB="GO_MF",Experiment="GTEX"),
  data.frame(newPE_PRvsNx_GTEXCC_max,DB="GO_CC",Experiment="GTEX"),
  data.frame(newPE_PRvsNx_GTEX_max,DB="All",Experiment="GTEX"),
  data.frame(newPE_PRvsNx_TCGABP_max,DB="GO_BP",Experiment="TCGA"),
  data.frame(newPE_PRvsNx_TCGAMF_max,DB="GO_MF",Experiment="TCGA"),
  data.frame(newPE_PRvsNx_TCGACC_max,DB="GO_CC",Experiment="TCGA"),
  data.frame(newPE_PRvsNx_TCGA_max,DB="All",Experiment="TCGA"))
newPE_PRvsNx_ggplot2$DB <- factor(newPE_PRvsNx_ggplot2$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))
p <- ggplot(newPE_PRvsNx_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")
pdf("newPE_PRxNx_histogram.pdf")
p
dev.off()
p <- ggplot(newPE_PRvsNx_ggplot2, aes(x=DB,y=Dist.max)) + geom_boxplot() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("newPE_PRxNx_boxplot.pdf")
p
dev.off()
p <- ggplot(newPE_PRvsNx_ggplot2, aes(x=DB,y=Dist.max)) + geom_violin() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("newPE_PRxNx_violin.pdf")
p
dev.off()

# wo MAX
newPE_PRvsNx_ggplot2woMax <- rbind(data.frame(newPE_PRvsNx_CCLEBP,DB="GO_BP",Experiment="CCLE"),
 data.frame(newPE_PRvsNx_CCLEMF,DB="GO_MF",Experiment="CCLE"),
 data.frame(newPE_PRvsNx_CCLECC,DB="GO_CC",Experiment="CCLE"),
 data.frame(newPE_PRvsNx_GTEXBP,DB="GO_BP",Experiment="GTEX"),
 data.frame(newPE_PRvsNx_GTEXMF,DB="GO_MF",Experiment="GTEX"),
 data.frame(newPE_PRvsNx_GTEXCC,DB="GO_CC",Experiment="GTEX"),
 data.frame(newPE_PRvsNx_TCGABP,DB="GO_BP",Experiment="TCGA"),
 data.frame(newPE_PRvsNx_TCGAMF,DB="GO_MF",Experiment="TCGA"),
 data.frame(newPE_PRvsNx_TCGACC,DB="GO_CC",Experiment="TCGA"))
newPE_PRvsNx_ggplot2woMax$DB <- factor(newPE_PRvsNx_ggplot2woMax$DB, levels=c("GO_BP","GO_MF","GO_CC"))
p <- ggplot(newPE_PRvsNx_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[c(1:3)]) + scale_fill_manual(values=brewer.pal(4,"Set1")[c(1:3)])
pdf("newPE_PRxNX_histogram_woMax.pdf")
p
dev.off()

##########################
# New known PE1 x CAFA
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_EN_CAFA/CCLE/",full.names=T)
newPE_PRvsCAFA_CCLEBP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsCAFA_CCLEBP <- tmp

  } else {

    newPE_PRvsCAFA_CCLEBP <- rbind(newPE_PRvsCAFA_CCLEBP, tmp)

  }

}
colnames(newPE_PRvsCAFA_CCLEBP) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsCAFA_CCLEBP$Dist <- as.numeric(paste(newPE_PRvsCAFA_CCLEBP$Dist))
newPE_PRvsCAFA_CCLEBP_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_CCLEBP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_EN_CAFA/CCLE/",full.names=T)
newPE_PRvsCAFA_CCLECC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsCAFA_CCLECC <- tmp

  } else {

    newPE_PRvsCAFA_CCLECC <- rbind(newPE_PRvsCAFA_CCLECC, tmp)

  }

}
colnames(newPE_PRvsCAFA_CCLECC) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsCAFA_CCLECC$Dist <- as.numeric(paste(newPE_PRvsCAFA_CCLECC$Dist))
newPE_PRvsCAFA_CCLECC_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_CCLECC[,c(1,5)], FUN=max, na.rm=T)

newPE_PRvsCAFA_CCLE <- rbind(newPE_PRvsCAFA_CCLEBP_max,newPE_PRvsCAFA_CCLECC_max)
newPE_PRvsCAFA_CCLE[,1] <- sapply(paste(newPE_PRvsCAFA_CCLE[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
newPE_PRvsCAFA_CCLE_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_CCLE, FUN=max, na.rm=T)
colnames(newPE_PRvsCAFA_CCLE_max)[2] <- "Dist.max"

# TCGA
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_EN_CAFA/TCGA/",full.names=T)
newPE_PRvsCAFA_TCGABP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsCAFA_TCGABP <- tmp

  } else {

    newPE_PRvsCAFA_TCGABP <- rbind(newPE_PRvsCAFA_TCGABP, tmp)

  }

}
colnames(newPE_PRvsCAFA_TCGABP) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsCAFA_TCGABP$Dist <- as.numeric(paste(newPE_PRvsCAFA_TCGABP$Dist))
newPE_PRvsCAFA_TCGABP_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_TCGABP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_EN_CAFA/TCGA/",full.names=T)
newPE_PRvsCAFA_TCGACC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsCAFA_TCGACC <- tmp

  } else {

    newPE_PRvsCAFA_TCGACC <- rbind(newPE_PRvsCAFA_TCGACC, tmp)

  }

}
colnames(newPE_PRvsCAFA_TCGACC) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsCAFA_TCGACC$Dist <- as.numeric(paste(newPE_PRvsCAFA_TCGACC$Dist))
newPE_PRvsCAFA_TCGACC_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_TCGACC[,c(1,5)], FUN=max, na.rm=T)

newPE_PRvsCAFA_TCGA <- rbind(newPE_PRvsCAFA_TCGABP_max,newPE_PRvsCAFA_TCGACC_max)
newPE_PRvsCAFA_TCGA[,1] <- sapply(paste(newPE_PRvsCAFA_TCGA[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
newPE_PRvsCAFA_TCGA_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_TCGA, FUN=max, na.rm=T)
colnames(newPE_PRvsCAFA_TCGA_max)[2] <- "Dist.max"

# GTEX
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_EN_CAFA/GTEX/",full.names=T)
newPE_PRvsCAFA_GTEXBP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsCAFA_GTEXBP <- tmp

  } else {

    newPE_PRvsCAFA_GTEXBP <- rbind(newPE_PRvsCAFA_GTEXBP, tmp)

  }

}
colnames(newPE_PRvsCAFA_GTEXBP) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsCAFA_GTEXBP$Dist <- as.numeric(paste(newPE_PRvsCAFA_GTEXBP$Dist))
newPE_PRvsCAFA_GTEXBP_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_GTEXBP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_EN_CAFA/GTEX/",full.names=T)
newPE_PRvsCAFA_GTEXCC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_PRvsCAFA_GTEXCC <- tmp

  } else {

    newPE_PRvsCAFA_GTEXCC <- rbind(newPE_PRvsCAFA_GTEXCC, tmp)

  }

}
colnames(newPE_PRvsCAFA_GTEXCC) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_PRvsCAFA_GTEXCC$Dist <- as.numeric(paste(newPE_PRvsCAFA_GTEXCC$Dist))
newPE_PRvsCAFA_GTEXCC_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_GTEXCC[,c(1,5)], FUN=max, na.rm=T)

newPE_PRvsCAFA_GTEX <- rbind(newPE_PRvsCAFA_GTEXBP_max,newPE_PRvsCAFA_GTEXCC_max)
newPE_PRvsCAFA_GTEX[,1] <- sapply(paste(newPE_PRvsCAFA_GTEX[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
newPE_PRvsCAFA_GTEX_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA_GTEX, FUN=max, na.rm=T)
colnames(newPE_PRvsCAFA_GTEX_max)[2] <- "Dist.max"

## ALL DBs and Datasets
newPE_PRvsCAFA <- rbind(newPE_PRvsCAFA_CCLE_max,newPE_PRvsCAFA_TCGA_max,newPE_PRvsCAFA_GTEX_max)
newPE_PRvsCAFA_max <- summaryBy(. ~ Prot, data=newPE_PRvsCAFA, FUN=max, na.rm=T)
colnames(newPE_PRvsCAFA_max)[2] <- "Dist.max"

newPE_PRvsCAFA_ggplot2 <- rbind(data.frame(newPE_PRvsCAFA_CCLEBP_max,DB="GO_BP",Experiment="CCLE"),
 data.frame(newPE_PRvsCAFA_CCLECC_max,DB="GO_CC",Experiment="CCLE"),
 data.frame(newPE_PRvsCAFA_CCLE_max,DB="All",Experiment="CCLE"),
 data.frame(newPE_PRvsCAFA_GTEXBP_max,DB="GO_BP",Experiment="GTEX"),
  data.frame(newPE_PRvsCAFA_GTEXCC_max,DB="GO_CC",Experiment="GTEX"),
  data.frame(newPE_PRvsCAFA_GTEX_max,DB="All",Experiment="GTEX"),
  data.frame(newPE_PRvsCAFA_TCGABP_max,DB="GO_BP",Experiment="TCGA"),
   data.frame(newPE_PRvsCAFA_TCGACC_max,DB="GO_CC",Experiment="TCGA"),
  data.frame(newPE_PRvsCAFA_TCGA_max,DB="All",Experiment="TCGA"))
newPE_PRvsCAFA_ggplot2$DB <- factor(newPE_PRvsCAFA_ggplot2$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))
p <- ggplot(newPE_PRvsCAFA_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[-2]) + scale_fill_manual(values=brewer.pal(4,"Set1")[-2])
pdf("newPE_PRxCAFA_histogram.pdf")
p
dev.off()
p <- ggplot(newPE_PRvsCAFA_ggplot2, aes(x=DB,y=Dist.max)) + geom_boxplot() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("newPE_PRxCAFA_boxplot.pdf")
p
dev.off()
p <- ggplot(newPE_PRvsCAFA_ggplot2, aes(x=DB,y=Dist.max)) + geom_violin() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("newPE_PRxCAFA_violin.pdf")
p
dev.off()

# wo MAX
newPE_PRvsCAFA_ggplot2woMax <- rbind(data.frame(newPE_PRvsCAFA_CCLEBP,DB="GO_BP",Experiment="CCLE"),
 data.frame(newPE_PRvsCAFA_CCLECC,DB="GO_CC",Experiment="CCLE"),
 data.frame(newPE_PRvsCAFA_GTEXBP,DB="GO_BP",Experiment="GTEX"),
  data.frame(newPE_PRvsCAFA_GTEXCC,DB="GO_CC",Experiment="GTEX"),
  data.frame(newPE_PRvsCAFA_TCGABP,DB="GO_BP",Experiment="TCGA"),
   data.frame(newPE_PRvsCAFA_TCGACC,DB="GO_CC",Experiment="TCGA"))
newPE_PRvsCAFA_ggplot2woMax$DB <- factor(newPE_PRvsCAFA_ggplot2woMax$DB, levels=c("GO_BP","GO_CC"))
p <- ggplot(newPE_PRvsCAFA_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[c(1,3)]) + scale_fill_manual(values=brewer.pal(4,"Set1")[c(1,3)])
pdf("newPE_PRxCAFA_histogram_woMax.pdf")
p
dev.off()

##########################
# New known PE1 NX x CAFA
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_CAFA_VS_NEXTPROT/",full.names=T)
newPE_NXvsCAFA_BP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_NXvsCAFA_BP <- tmp

  } else {

    newPE_NXvsCAFA_BP <- rbind(newPE_NXvsCAFA_BP, tmp)

  }

}
colnames(newPE_NXvsCAFA_BP) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_NXvsCAFA_BP$Dist <- as.numeric(paste(newPE_NXvsCAFA_BP$Dist))
newPE_NXvsCAFA_BP_max <- summaryBy(. ~ Prot, data=newPE_NXvsCAFA_BP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/uPES_PE1_known_ENRIQUECIDAS_PR_CAFA_VS_NEXTPROT/",full.names=T)
newPE_NXvsCAFA_CC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    newPE_NXvsCAFA_CC <- tmp

  } else {

    newPE_NXvsCAFA_CC <- rbind(newPE_NXvsCAFA_CC, tmp)

  }

}
colnames(newPE_NXvsCAFA_CC) <- c("Prot","GO1","GO2","Ont","Dist")
newPE_NXvsCAFA_CC$Dist <- as.numeric(paste(newPE_NXvsCAFA_CC$Dist))
newPE_NXvsCAFA_CC_max <- summaryBy(. ~ Prot, data=newPE_NXvsCAFA_CC[,c(1,5)], FUN=max, na.rm=T)

newPE_NXvsCAFA <- rbind(newPE_NXvsCAFA_BP_max,newPE_NXvsCAFA_CC_max)
newPE_NXvsCAFA[,1] <- sapply(paste(newPE_NXvsCAFA[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
newPE_NXvsCAFA_max <- summaryBy(. ~ Prot, data=newPE_NXvsCAFA, FUN=max, na.rm=T)
colnames(newPE_NXvsCAFA_max)[2] <- "Dist.max"

newPE_NXvsCAFA_ggplot2 <- rbind(data.frame(newPE_NXvsCAFA_BP_max,DB="GO_BP",Experiment="CAFAxNX"),
 data.frame(newPE_NXvsCAFA_CC_max,DB="GO_CC",Experiment="CAFAxNX"),
 data.frame(newPE_NXvsCAFA_max,DB="All",Experiment="CAFAxNX"))
newPE_NXvsCAFA_ggplot2$DB <- factor(newPE_NXvsCAFA_ggplot2$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(newPE_NXvsCAFA_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(. ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[-2]) + scale_fill_manual(values=brewer.pal(4,"Set1")[-2])
pdf("newPE_CAFAxNX_histogram.pdf")
p
dev.off()
p <- ggplot(newPE_NXvsCAFA_ggplot2, aes(x=DB,y=Dist.max)) + geom_boxplot() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank())
pdf("newPE_CAFAxNX_boxplot.pdf")
p
dev.off()
p <- ggplot(newPE_NXvsCAFA_ggplot2, aes(x=DB,y=Dist.max)) + geom_violin() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank())
pdf("newPE_CAFAxNX_violin.pdf")
p
dev.off()

# wo MAX
newPE_NXvsCAFA_ggplot2woMax <- rbind(data.frame(newPE_NXvsCAFA_BP_max,DB="GO_BP",Experiment="CAFAxNX"),
  data.frame(newPE_NXvsCAFA_CC_max,DB="GO_CC",Experiment="CAFAxNX"),
 data.frame(newPE_NXvsCAFA_max,DB="All",Experiment="CAFAxNX"))
newPE_NXvsCAFA_ggplot2woMax$DB <- factor(newPE_NXvsCAFA_ggplot2woMax$DB, levels=c("GO_BP","GO_CC","All"))
p <- ggplot(newPE_NXvsCAFA_ggplot2woMax, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[c(1,3,4)]) + scale_fill_manual(values=brewer.pal(4,"Set1")[c(1,3,4)])
pdf("newPE_CAFAxNX_histogram_woMax.pdf")
p
dev.off()


##########################
# 10% PE1 
##########################

##########################
# 10% PE1 x nextProt
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE/",full.names=T)
PE10_PRvsNx_CCLEBP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_CCLEBP <- tmp

  } else {

    PE10_PRvsNx_CCLEBP <- rbind(PE10_PRvsNx_CCLEBP, tmp)

  }

}
colnames(PE10_PRvsNx_CCLEBP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_CCLEBP$Dist <- as.numeric(paste(PE10_PRvsNx_CCLEBP$Dist))
PE10_PRvsNx_CCLEBP_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_CCLEBP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE/",full.names=T)
PE10_PRvsNx_CCLEMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_CCLEMF <- tmp

  } else {

    PE10_PRvsNx_CCLEMF <- rbind(PE10_PRvsNx_CCLEMF, tmp)

  }

}
colnames(PE10_PRvsNx_CCLEMF) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_CCLEMF$Dist <- as.numeric(paste(PE10_PRvsNx_CCLEMF$Dist))
PE10_PRvsNx_CCLEMF_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_CCLEMF[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsNx_CCLEMF_max[PE10_PRvsNx_CCLEMF_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE/",full.names=T)
PE10_PRvsNx_CCLECC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_CCLECC <- tmp

  } else {

    PE10_PRvsNx_CCLECC <- rbind(PE10_PRvsNx_CCLECC, tmp)

  }

}
colnames(PE10_PRvsNx_CCLECC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_CCLECC$Dist <- as.numeric(paste(PE10_PRvsNx_CCLECC$Dist))
PE10_PRvsNx_CCLECC_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_CCLECC[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsNx_CCLECC_max[PE10_PRvsNx_CCLECC_max==-Inf] <- NA

PE10_PRvsNx_CCLE <- rbind(PE10_PRvsNx_CCLEBP_max,PE10_PRvsNx_CCLEMF_max,PE10_PRvsNx_CCLECC_max)
PE10_PRvsNx_CCLE[,1] <- sapply(paste(PE10_PRvsNx_CCLE[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsNx_CCLE_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_CCLE, FUN=max, na.rm=T)
colnames(PE10_PRvsNx_CCLE_max)[2] <- "Dist.max"

# TCGA
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA/",full.names=T)
PE10_PRvsNx_TCGABP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_TCGABP <- tmp

  } else {

    PE10_PRvsNx_TCGABP <- rbind(PE10_PRvsNx_TCGABP, tmp)

  }

}
colnames(PE10_PRvsNx_TCGABP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_TCGABP$Dist <- as.numeric(paste(PE10_PRvsNx_TCGABP$Dist))
PE10_PRvsNx_TCGABP_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_TCGABP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA/",full.names=T)
PE10_PRvsNx_TCGAMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_TCGAMF <- tmp

  } else {

    PE10_PRvsNx_TCGAMF <- rbind(PE10_PRvsNx_TCGAMF, tmp)

  }

}
colnames(PE10_PRvsNx_TCGAMF) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_TCGAMF$Dist <- as.numeric(paste(PE10_PRvsNx_TCGAMF$Dist))
PE10_PRvsNx_TCGAMF_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_TCGAMF[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA/",full.names=T)
PE10_PRvsNx_TCGACC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_TCGACC <- tmp

  } else {

    PE10_PRvsNx_TCGACC <- rbind(PE10_PRvsNx_TCGACC, tmp)

  }

}
colnames(PE10_PRvsNx_TCGACC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_TCGACC$Dist <- as.numeric(paste(PE10_PRvsNx_TCGACC$Dist))
PE10_PRvsNx_TCGACC_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_TCGACC[,c(1,5)], FUN=max, na.rm=T)

PE10_PRvsNx_TCGA <- rbind(PE10_PRvsNx_TCGABP_max,PE10_PRvsNx_TCGAMF_max,PE10_PRvsNx_TCGACC_max)
PE10_PRvsNx_TCGA[,1] <- sapply(paste(PE10_PRvsNx_TCGA[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsNx_TCGA_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_TCGA, FUN=max, na.rm=T)
colnames(PE10_PRvsNx_TCGA_max)[2] <- "Dist.max"

# GTEX
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX/",full.names=T)
PE10_PRvsNx_GTEXBP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_GTEXBP <- tmp

  } else {

    PE10_PRvsNx_GTEXBP <- rbind(PE10_PRvsNx_GTEXBP, tmp)

  }

}
colnames(PE10_PRvsNx_GTEXBP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_GTEXBP$Dist <- as.numeric(paste(PE10_PRvsNx_GTEXBP$Dist))
PE10_PRvsNx_GTEXBP_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_GTEXBP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX/",full.names=T)
PE10_PRvsNx_GTEXMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_GTEXMF <- tmp

  } else {

    PE10_PRvsNx_GTEXMF <- rbind(PE10_PRvsNx_GTEXMF, tmp)

  }

}
colnames(PE10_PRvsNx_GTEXMF) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_GTEXMF$Dist <- as.numeric(paste(PE10_PRvsNx_GTEXMF$Dist))
PE10_PRvsNx_GTEXMF_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_GTEXMF[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX/",full.names=T)
PE10_PRvsNx_GTEXCC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_GTEXCC <- tmp

  } else {

    PE10_PRvsNx_GTEXCC <- rbind(PE10_PRvsNx_GTEXCC, tmp)

  }

}
colnames(PE10_PRvsNx_GTEXCC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_GTEXCC$Dist <- as.numeric(paste(PE10_PRvsNx_GTEXCC$Dist))
PE10_PRvsNx_GTEXCC_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_GTEXCC[,c(1,5)], FUN=max, na.rm=T)

PE10_PRvsNx_GTEX <- rbind(PE10_PRvsNx_GTEXBP_max,PE10_PRvsNx_GTEXMF_max,PE10_PRvsNx_GTEXCC_max)
PE10_PRvsNx_GTEX[,1] <- sapply(paste(PE10_PRvsNx_GTEX[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsNx_GTEX_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_GTEX, FUN=max, na.rm=T)
colnames(PE10_PRvsNx_GTEX_max)[2] <- "Dist.max"

# wo MAX
PE10_PRvsNX_ggplot2woMax <- rbind(data.frame(PE10_PRvsNx_CCLEBP,DB="GO_BP",Experiment="CCLE"),
  data.frame(PE10_PRvsNx_CCLEMF,DB="GO_MF",Experiment="CCLE"),
 data.frame(PE10_PRvsNx_CCLECC,DB="GO_CC",Experiment="CCLE"),
  data.frame(PE10_PRvsNx_GTEXBP,DB="GO_BP",Experiment="GTEX"),
 data.frame(PE10_PRvsNx_GTEXMF,DB="GO_MF",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_GTEXCC,DB="GO_CC",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_TCGABP,DB="GO_BP",Experiment="TCGA"),
  data.frame(PE10_PRvsNx_TCGAMF,DB="GO_MF",Experiment="TCGA"),
   data.frame(PE10_PRvsNx_TCGACC,DB="GO_CC",Experiment="TCGA"))
PE10_PRvsNX_ggplot2woMax$DB <- factor(PE10_PRvsNX_ggplot2woMax$DB, levels=c("GO_BP","GO_MF","GO_CC"))

p <- ggplot(PE10_PRvsNX_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[1:3]) + scale_fill_manual(values=brewer.pal(4,"Set1")[1:3])
pdf("PE10_PRxNX_histogram_woMax.pdf")
p
dev.off()



## ALL DBs and Datasets
PE10_PRvsNx <- rbind(PE10_PRvsNx_CCLE_max,PE10_PRvsNx_TCGA_max,PE10_PRvsNx_GTEX_max)
PE10_PRvsNx_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx, FUN=max, na.rm=T)
colnames(PE10_PRvsNx_max)[2] <- "Dist.max"

PE10_PRvsNx_ggplot2 <- rbind(data.frame(PE10_PRvsNx_CCLEBP_max,DB="GO_BP",Experiment="CCLE"),
  data.frame(PE10_PRvsNx_CCLEMF_max,DB="GO_MF",Experiment="CCLE"),
 data.frame(PE10_PRvsNx_CCLECC_max,DB="GO_CC",Experiment="CCLE"),
 data.frame(PE10_PRvsNx_CCLE_max,DB="All",Experiment="CCLE"),
 data.frame(PE10_PRvsNx_GTEXBP_max,DB="GO_BP",Experiment="GTEX"),
 data.frame(PE10_PRvsNx_GTEXMF_max,DB="GO_MF",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_GTEXCC_max,DB="GO_CC",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_GTEX_max,DB="All",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_TCGABP_max,DB="GO_BP",Experiment="TCGA"),
  data.frame(PE10_PRvsNx_TCGAMF_max,DB="GO_MF",Experiment="TCGA"),
  data.frame(PE10_PRvsNx_TCGACC_max,DB="GO_CC",Experiment="TCGA"),
  data.frame(PE10_PRvsNx_TCGA_max,DB="All",Experiment="TCGA"))
PE10_PRvsNx_ggplot2$DB <- factor(PE10_PRvsNx_ggplot2$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(PE10_PRvsNx_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")
pdf("PE10_PRxNx_histogram.pdf")
p
dev.off()
p <- ggplot(PE10_PRvsNx_ggplot2, aes(x=DB,y=Dist.max)) + geom_boxplot() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PE10_PRxNx_boxplot.pdf")
p
dev.off()
p <- ggplot(PE10_PRvsNx_ggplot2, aes(x=DB,y=Dist.max)) + geom_violin() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PE10_PRxNx_violin.pdf")
p
dev.off()


##########################
# 10% PE1 x CAFA
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE/",full.names=T)
PE10_PRvsCAFA_CCLEBP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_CCLEBP <- tmp

  } else {

    PE10_PRvsCAFA_CCLEBP <- rbind(PE10_PRvsCAFA_CCLEBP, tmp)

  }

}
colnames(PE10_PRvsCAFA_CCLEBP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_CCLEBP$Dist <- as.numeric(paste(PE10_PRvsCAFA_CCLEBP$Dist))
PE10_PRvsCAFA_CCLEBP_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLEBP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE/",full.names=T)
PE10_PRvsCAFA_CCLEMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_CCLEMF <- tmp

  } else {

    PE10_PRvsCAFA_CCLEMF <- rbind(PE10_PRvsCAFA_CCLEMF, tmp)

  }

}
colnames(PE10_PRvsCAFA_CCLEMF) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_CCLEMF$Dist <- as.numeric(paste(PE10_PRvsCAFA_CCLEMF$Dist))
PE10_PRvsCAFA_CCLEMF_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLEMF[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE/",full.names=T)
PE10_PRvsCAFA_CCLECC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_CCLECC <- tmp

  } else {

    PE10_PRvsCAFA_CCLECC <- rbind(PE10_PRvsCAFA_CCLECC, tmp)

  }

}
colnames(PE10_PRvsCAFA_CCLECC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_CCLECC$Dist <- as.numeric(paste(PE10_PRvsCAFA_CCLECC$Dist))
PE10_PRvsCAFA_CCLECC_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLECC[,c(1,5)], FUN=max, na.rm=T)

PE10_PRvsCAFA_CCLE <- rbind(PE10_PRvsCAFA_CCLEBP_max, PE10_PRvsCAFA_CCLEMF_max, PE10_PRvsCAFA_CCLECC_max)
PE10_PRvsCAFA_CCLE[,1] <- sapply(paste(PE10_PRvsCAFA_CCLE[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_CCLE_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLE, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFA_CCLE_max)[2] <- "Dist.max"

# TCGA
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA/",full.names=T)
PE10_PRvsCAFA_TCGABP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_TCGABP <- tmp

  } else {

    PE10_PRvsCAFA_TCGABP <- rbind(PE10_PRvsCAFA_TCGABP, tmp)

  }

}
colnames(PE10_PRvsCAFA_TCGABP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_TCGABP$Dist <- as.numeric(paste(PE10_PRvsCAFA_TCGABP$Dist))
PE10_PRvsCAFA_TCGABP_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGABP[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_TCGABP_max[PE10_PRvsCAFA_TCGABP_max==-Inf] <- NA

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA/",full.names=T)
PE10_PRvsCAFA_TCGAMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_TCGAMF <- tmp

  } else {

    PE10_PRvsCAFA_TCGAMF <- rbind(PE10_PRvsCAFA_TCGAMF, tmp)

  }

}
colnames(PE10_PRvsCAFA_TCGAMF) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_TCGAMF$Dist <- as.numeric(paste(PE10_PRvsCAFA_TCGAMF$Dist))
PE10_PRvsCAFA_TCGAMF_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGAMF[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_TCGAMF_max[PE10_PRvsCAFA_TCGAMF_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA/",full.names=T)
PE10_PRvsCAFA_TCGACC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_TCGACC <- tmp

  } else {

    PE10_PRvsCAFA_TCGACC <- rbind(PE10_PRvsCAFA_TCGACC, tmp)

  }

}
colnames(PE10_PRvsCAFA_TCGACC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_TCGACC$Dist <- as.numeric(paste(PE10_PRvsCAFA_TCGACC$Dist))
PE10_PRvsCAFA_TCGACC_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGACC[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_TCGACC_max[PE10_PRvsCAFA_TCGACC_max==-Inf] <- NA

PE10_PRvsCAFA_TCGA <- rbind(PE10_PRvsCAFA_TCGABP_max,PE10_PRvsCAFA_TCGAMF_max,PE10_PRvsCAFA_TCGACC_max)
PE10_PRvsCAFA_TCGA[,1] <- sapply(paste(PE10_PRvsCAFA_TCGA[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_TCGA_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGA, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFA_TCGA_max)[2] <- "Dist.max"

# GTEX
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX/",full.names=T)
PE10_PRvsCAFA_GTEXBP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_GTEXBP <- tmp

  } else {

    PE10_PRvsCAFA_GTEXBP <- rbind(PE10_PRvsCAFA_GTEXBP, tmp)

  }

}
colnames(PE10_PRvsCAFA_GTEXBP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_GTEXBP$Dist <- as.numeric(paste(PE10_PRvsCAFA_GTEXBP$Dist))
PE10_PRvsCAFA_GTEXBP_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEXBP[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_GTEXBP_max[PE10_PRvsCAFA_GTEXBP_max==-Inf] <- NA

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX/",full.names=T)
PE10_PRvsCAFA_GTEXMF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_GTEXMF <- tmp

  } else {

    PE10_PRvsCAFA_GTEXMF <- rbind(PE10_PRvsCAFA_GTEXMF, tmp)

  }
}
colnames(PE10_PRvsCAFA_GTEXMF) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_GTEXMF$Dist <- as.numeric(paste(PE10_PRvsCAFA_GTEXMF$Dist))
PE10_PRvsCAFA_GTEXMF_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEXMF[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_GTEXMF_max[PE10_PRvsCAFA_GTEXMF_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX/",full.names=T)
PE10_PRvsCAFA_GTEXCC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_GTEXCC <- tmp

  } else {

    PE10_PRvsCAFA_GTEXCC <- rbind(PE10_PRvsCAFA_GTEXCC, tmp)

  }

}
colnames(PE10_PRvsCAFA_GTEXCC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_GTEXCC$Dist <- as.numeric(paste(PE10_PRvsCAFA_GTEXCC$Dist))
PE10_PRvsCAFA_GTEXCC_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEXCC[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_GTEXCC_max[PE10_PRvsCAFA_GTEXCC_max==-Inf] <- NA

PE10_PRvsCAFA_GTEX <- rbind(PE10_PRvsCAFA_GTEXBP_max,PE10_PRvsCAFA_GTEXMF_max,PE10_PRvsCAFA_GTEXCC_max)
PE10_PRvsCAFA_GTEX[,1] <- sapply(paste(PE10_PRvsCAFA_GTEX[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_GTEX_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEX, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFA_GTEX_max)[2] <- "Dist.max"
PE10_PRvsCAFA_GTEX_max[PE10_PRvsCAFA_GTEX_max==-Inf] <- NA

## ALL DBs and Datasets
PE10_PRvsCAFA <- rbind(PE10_PRvsCAFA_CCLE_max,PE10_PRvsCAFA_TCGA_max,PE10_PRvsCAFA_GTEX_max)
PE10_PRvsCAFA_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFA_max)[2] <- "Dist.max"

PE10_PRvsCAFA_ggplot2 <- rbind(data.frame(PE10_PRvsCAFA_CCLEBP_max,DB="GO_BP",Experiment="CCLE"),
  data.frame(PE10_PRvsCAFA_CCLEMF_max,DB="GO_MF",Experiment="CCLE"),
 data.frame(PE10_PRvsCAFA_CCLECC_max,DB="GO_CC",Experiment="CCLE"),
 data.frame(PE10_PRvsCAFA_CCLE_max,DB="All",Experiment="CCLE"),
 data.frame(PE10_PRvsCAFA_GTEXBP_max,DB="GO_BP",Experiment="GTEX"),
 data.frame(PE10_PRvsCAFA_GTEXMF_max,DB="GO_MF",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_GTEXCC_max,DB="GO_CC",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_GTEX_max,DB="All",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_TCGABP_max,DB="GO_BP",Experiment="TCGA"),
  data.frame(PE10_PRvsCAFA_TCGAMF_max,DB="GO_MF",Experiment="TCGA"),
   data.frame(PE10_PRvsCAFA_TCGACC_max,DB="GO_CC",Experiment="TCGA"),
  data.frame(PE10_PRvsCAFA_TCGA_max,DB="All",Experiment="TCGA"))
PE10_PRvsCAFA_ggplot2$DB <- factor(PE10_PRvsCAFA_ggplot2$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(PE10_PRvsCAFA_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")) + scale_fill_manual(values=brewer.pal(4,"Set1"))
pdf("PE10_PRxCAFA_histogram.pdf")
p
dev.off()
p <- ggplot(PE10_PRvsCAFA_ggplot2, aes(x=DB,y=Dist.max)) + geom_boxplot() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PE10_PRxCAFA_boxplot.pdf")
p
dev.off()
p <- ggplot(PE10_PRvsCAFA_ggplot2, aes(x=DB,y=Dist.max)) + geom_violin() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PE10_PRxCAFA_violin.pdf")
p
dev.off()

# wo MAX
PE10_PRvsCAFA_ggplot2woMax <- rbind(data.frame(PE10_PRvsCAFA_CCLEBP,DB="GO_BP",Experiment="CCLE"),
  data.frame(PE10_PRvsCAFA_CCLEMF,DB="GO_MF",Experiment="CCLE"),
 data.frame(PE10_PRvsCAFA_CCLECC,DB="GO_CC",Experiment="CCLE"),
  data.frame(PE10_PRvsCAFA_GTEXBP,DB="GO_BP",Experiment="GTEX"),
 data.frame(PE10_PRvsCAFA_GTEXMF,DB="GO_MF",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_GTEXCC,DB="GO_CC",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_TCGABP,DB="GO_BP",Experiment="TCGA"),
  data.frame(PE10_PRvsCAFA_TCGAMF,DB="GO_MF",Experiment="TCGA"),
   data.frame(PE10_PRvsCAFA_TCGACC,DB="GO_CC",Experiment="TCGA"))
PE10_PRvsCAFA_ggplot2woMax$DB <- factor(PE10_PRvsCAFA_ggplot2woMax$DB, levels=c("GO_BP","GO_MF","GO_CC"))

p <- ggplot(PE10_PRvsCAFA_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[1:3]) + scale_fill_manual(values=brewer.pal(4,"Set1")[1:3])
pdf("PE10_PRxCAFA_histogram_woMax.pdf")
p
dev.off()

##########################
# PE1 10% NX x CAFA
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_CAFA_vs_NEXTPROT/",full.names=T)
PE10_NXvsCAFA_BP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_NXvsCAFA_BP <- tmp

  } else {

    PE10_NXvsCAFA_BP <- rbind(PE10_NXvsCAFA_BP, tmp)

  }

}
colnames(PE10_NXvsCAFA_BP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_NXvsCAFA_BP$Dist <- as.numeric(paste(PE10_NXvsCAFA_BP$Dist))
PE10_NXvsCAFA_BP_max <- summaryBy(. ~ Prot, data=PE10_NXvsCAFA_BP[,c(1,5)], FUN=max, na.rm=T)
PE10_NXvsCAFA_BP_max[PE10_NXvsCAFA_BP_max==-Inf] <- NA 

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_CAFA_vs_NEXTPROT/",full.names=T)
PE10_NXvsCAFA_MF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_NXvsCAFA_MF <- tmp

  } else {

    PE10_NXvsCAFA_MF <- rbind(PE10_NXvsCAFA_MF, tmp)

  }

}
colnames(PE10_NXvsCAFA_MF) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_NXvsCAFA_MF$Dist <- as.numeric(paste(PE10_NXvsCAFA_MF$Dist))
PE10_NXvsCAFA_MF_max <- summaryBy(. ~ Prot, data=PE10_NXvsCAFA_MF[,c(1,5)], FUN=max, na.rm=T)
PE10_NXvsCAFA_MF_max[PE10_NXvsCAFA_MF_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_CAFA_vs_NEXTPROT/",full.names=T)
PE10_NXvsCAFA_CC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_NXvsCAFA_CC <- tmp

  } else {

    PE10_NXvsCAFA_CC <- rbind(PE10_NXvsCAFA_CC, tmp)

  }

}
colnames(PE10_NXvsCAFA_CC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_NXvsCAFA_CC$Dist <- as.numeric(paste(PE10_NXvsCAFA_CC$Dist))
PE10_NXvsCAFA_CC_max <- summaryBy(. ~ Prot, data=PE10_NXvsCAFA_CC[,c(1,5)], FUN=max, na.rm=T)
PE10_NXvsCAFA_CC_max[PE10_NXvsCAFA_CC_max==-Inf] <- NA

PE10_NXvsCAFA <- rbind(PE10_NXvsCAFA_BP_max,PE10_NXvsCAFA_MF_max,PE10_NXvsCAFA_CC_max)
PE10_NXvsCAFA[,1] <- sapply(paste(PE10_NXvsCAFA[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_NXvsCAFA_max <- summaryBy(. ~ Prot, data=PE10_NXvsCAFA, FUN=max, na.rm=T)
colnames(PE10_NXvsCAFA_max)[2] <- "Dist.max"

PE10_NXvsCAFA_ggplot2 <- rbind(data.frame(PE10_NXvsCAFA_BP_max,DB="GO_BP",Experiment="CAFAxNX"),
  data.frame(PE10_NXvsCAFA_MF_max,DB="GO_MF",Experiment="CAFAxNX"),
 data.frame(PE10_NXvsCAFA_CC_max,DB="GO_CC",Experiment="CAFAxNX"),
 data.frame(PE10_NXvsCAFA_max,DB="All",Experiment="CAFAxNX"))
PE10_NXvsCAFA_ggplot2$DB <- factor(PE10_NXvsCAFA_ggplot2$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(PE10_NXvsCAFA_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(. ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")) + scale_fill_manual(values=brewer.pal(4,"Set1"))
pdf("PE10_CAFAxNX_histogram.pdf")
p
dev.off()
p <- ggplot(PE10_NXvsCAFA_ggplot2, aes(x=DB,y=Dist.max)) + geom_boxplot() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank())
pdf("PE10_CAFAxNX_boxplot.pdf")
p
dev.off()
p <- ggplot(PE10_NXvsCAFA_ggplot2, aes(x=DB,y=Dist.max)) + geom_violin() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank())
pdf("PE10_CAFAxNX_violin.pdf")
p
dev.off()

# wo MAX
PE10_NXvsCAFA_ggplot2woMax <- rbind(data.frame(PE10_NXvsCAFA_BP,DB="GO_BP",Experiment="CAFAxNX"),
  data.frame(PE10_NXvsCAFA_MF,DB="GO_MF",Experiment="CAFAxNX"),
 data.frame(PE10_NXvsCAFA_CC,DB="GO_CC",Experiment="CAFAxNX"))
PE10_NXvsCAFA_ggplot2woMax$DB <- factor(PE10_NXvsCAFA_ggplot2woMax$DB, levels=c("GO_BP","GO_MF","GO_CC"))

p <- ggplot(PE10_NXvsCAFA_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(. ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[1:3]) + scale_fill_manual(values=brewer.pal(4,"Set1")[1:3])
pdf("PE10_CAFAxNX_histogram_woMax.pdf")
p
dev.off()


benchmarking_ggplot2 <- rbind(data.frame(newPE_NXvsCAFA_max,Benchmark="NXvsCAFA", ProteinSet="newPE1"),data.frame(newPE_PRvsNx_max,Benchmark="PRvsNX", ProteinSet="newPE1"), data.frame(newPE_PRvsCAFA_max,Benchmark="PRvsCAFA", ProteinSet="newPE1"),data.frame(PE10_NXvsCAFA_max,Benchmark="NXvsCAFA", ProteinSet="PE1_10%"),data.frame(PE10_PRvsNx_max,Benchmark="PRvsNX", ProteinSet="PE1_10%"), data.frame(PE10_PRvsCAFA_max,Benchmark="PRvsCAFA", ProteinSet="PE1_10%"))
p <- ggplot(benchmarking_ggplot2, aes(x=Dist.max, col=Benchmark, fill=Benchmark)) + geom_histogram() + facet_grid(ProteinSet ~ Benchmark, scales="free") + theme_bw()
pdf("Benchmarking_allDBallExp_histogram.pdf")
p
dev.off()



##########################
# 10% PE1 RANDOM
##########################

##########################
# 10% PE1 x nextProt
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_CCLEBPr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_CCLEBPr <- tmp

  } else {

    PE10_PRvsNx_CCLEBPr <- rbind(PE10_PRvsNx_CCLEBPr, tmp)

  }

}
colnames(PE10_PRvsNx_CCLEBPr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_CCLEBPr$Dist <- as.numeric(paste(PE10_PRvsNx_CCLEBPr$Dist))
PE10_PRvsNx_CCLEBPr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_CCLEBPr[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_CCLEMFr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_CCLEMFr <- tmp

  } else {

    PE10_PRvsNx_CCLEMFr <- rbind(PE10_PRvsNx_CCLEMFr, tmp)

  }

}
colnames(PE10_PRvsNx_CCLEMFr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_CCLEMFr$Dist <- as.numeric(paste(PE10_PRvsNx_CCLEMFr$Dist))
PE10_PRvsNx_CCLEMFr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_CCLEMFr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsNx_CCLEMFr_max[PE10_PRvsNx_CCLEMFr_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_CCLECCr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_CCLECCr <- tmp

  } else {

    PE10_PRvsNx_CCLECCr <- rbind(PE10_PRvsNx_CCLECCr, tmp)

  }
}
colnames(PE10_PRvsNx_CCLECCr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_CCLECCr$Dist <- as.numeric(paste(PE10_PRvsNx_CCLECCr$Dist))
PE10_PRvsNx_CCLECCr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_CCLECCr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsNx_CCLECCr_max[PE10_PRvsNx_CCLECCr_max==-Inf] <- NA

PE10_PRvsNx_CCLEr <- rbind(PE10_PRvsNx_CCLEBPr_max,PE10_PRvsNx_CCLEMFr_max,PE10_PRvsNx_CCLECCr_max)
PE10_PRvsNx_CCLEr[,1] <- sapply(paste(PE10_PRvsNx_CCLEr[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsNx_CCLEr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_CCLEr, FUN=max, na.rm=T)
colnames(PE10_PRvsNx_CCLEr_max)[2] <- "Dist.max"

# TCGA
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_TCGABPr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_TCGABPr <- tmp

  } else {

    PE10_PRvsNx_TCGABPr <- rbind(PE10_PRvsNx_TCGABPr, tmp)

  }

}
colnames(PE10_PRvsNx_TCGABPr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_TCGABPr$Dist <- as.numeric(paste(PE10_PRvsNx_TCGABPr$Dist))
PE10_PRvsNx_TCGABPr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_TCGABPr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsNx_TCGABPr_max[PE10_PRvsNx_TCGABPr_max==-Inf] <- NA

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_TCGAMFr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_TCGAMFr <- tmp

  } else {

    PE10_PRvsNx_TCGAMFr <- rbind(PE10_PRvsNx_TCGAMFr, tmp)

  }

}
colnames(PE10_PRvsNx_TCGAMFr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_TCGAMFr$Dist <- as.numeric(paste(PE10_PRvsNx_TCGAMFr$Dist))
PE10_PRvsNx_TCGAMFr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_TCGAMFr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsNx_TCGAMFr_max[PE10_PRvsNx_TCGAMFr_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_TCGACCr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_TCGACCr <- tmp

  } else {

    PE10_PRvsNx_TCGACCr <- rbind(PE10_PRvsNx_TCGACCr, tmp)

  }

}
colnames(PE10_PRvsNx_TCGACCr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_TCGACCr$Dist <- as.numeric(paste(PE10_PRvsNx_TCGACCr$Dist))
PE10_PRvsNx_TCGACCr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_TCGACCr[,c(1,5)], FUN=max, na.rm=T)

PE10_PRvsNx_TCGAr <- rbind(PE10_PRvsNx_TCGABPr_max,PE10_PRvsNx_TCGAMFr_max,PE10_PRvsNx_TCGACCr_max)
PE10_PRvsNx_TCGAr[,1] <- sapply(paste(PE10_PRvsNx_TCGAr[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsNx_TCGAr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_TCGAr, FUN=max, na.rm=T)
colnames(PE10_PRvsNx_TCGAr_max)[2] <- "Dist.max"

# GTEX
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_GTEXBPr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_GTEXBPr <- tmp

  } else {

    PE10_PRvsNx_GTEXBPr <- rbind(PE10_PRvsNx_GTEXBPr, tmp)

  }

}
colnames(PE10_PRvsNx_GTEXBPr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_GTEXBPr$Dist <- as.numeric(paste(PE10_PRvsNx_GTEXBPr$Dist))
PE10_PRvsNx_GTEXBPr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_GTEXBPr[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_GTEXMFr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_GTEXMFr <- tmp

  } else {

    PE10_PRvsNx_GTEXMFr <- rbind(PE10_PRvsNx_GTEXMFr, tmp)

  }

}
colnames(PE10_PRvsNx_GTEXMFr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_GTEXMFr$Dist <- as.numeric(paste(PE10_PRvsNx_GTEXMFr$Dist))
PE10_PRvsNx_GTEXMFr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_GTEXMFr[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsNx_GTEXCCr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_GTEXCCr <- tmp

  } else {

    PE10_PRvsNx_GTEXCCr <- rbind(PE10_PRvsNx_GTEXCCr, tmp)

  }

}
colnames(PE10_PRvsNx_GTEXCCr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_GTEXCCr$Dist <- as.numeric(paste(PE10_PRvsNx_GTEXCCr$Dist))
PE10_PRvsNx_GTEXCCr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_GTEXCCr[,c(1,5)], FUN=max, na.rm=T)

PE10_PRvsNx_GTEXr <- rbind(PE10_PRvsNx_GTEXBPr_max,PE10_PRvsNx_GTEXMFr_max,PE10_PRvsNx_GTEXCCr_max)
PE10_PRvsNx_GTEXr[,1] <- sapply(paste(PE10_PRvsNx_GTEXr[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsNx_GTEXr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_GTEXr, FUN=max, na.rm=T)
colnames(PE10_PRvsNx_GTEXr_max)[2] <- "Dist.max"

## ALL DBs and Datasets
PE10_PRvsNxr <- rbind(PE10_PRvsNx_CCLEr_max,PE10_PRvsNx_TCGAr_max,PE10_PRvsNx_GTEXr_max)
PE10_PRvsNxr_max <- summaryBy(. ~ Prot, data=PE10_PRvsNxr, FUN=max, na.rm=T)
colnames(PE10_PRvsNxr_max)[2] <- "Dist.max"

PE10_PRvsNx_ggplot2_r <- rbind(data.frame(PE10_PRvsNx_CCLEBPr_max,DB="GO_BP",Experiment="CCLE"),
  data.frame(PE10_PRvsNx_CCLEMFr_max,DB="GO_MF",Experiment="CCLE"),
 data.frame(PE10_PRvsNx_CCLECCr_max,DB="GO_CC",Experiment="CCLE"),
 data.frame(PE10_PRvsNx_CCLEr_max,DB="All",Experiment="CCLE"),
 data.frame(PE10_PRvsNx_GTEXBPr_max,DB="GO_BP",Experiment="GTEX"),
 data.frame(PE10_PRvsNx_GTEXMFr_max,DB="GO_MF",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_GTEXCCr_max,DB="GO_CC",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_GTEXr_max,DB="All",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_TCGABPr_max,DB="GO_BP",Experiment="TCGA"),
  data.frame(PE10_PRvsNx_TCGAMFr_max,DB="GO_MF",Experiment="TCGA"),
  data.frame(PE10_PRvsNx_TCGACCr_max,DB="GO_CC",Experiment="TCGA"),
  data.frame(PE10_PRvsNx_TCGAr_max,DB="All",Experiment="TCGA"))
PE10_PRvsNx_ggplot2_r$DB <- factor(PE10_PRvsNx_ggplot2_r$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(PE10_PRvsNx_ggplot2_r, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")
pdf("PE10_PRxNx_Random2iteration1_histogram.pdf")
p
dev.off()
p <- ggplot(PE10_PRvsNx_ggplot2_r, aes(x=DB,y=Dist.max)) + geom_boxplot() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PE10_PRxNx_Random2_boxplot.pdf")
p
dev.off()
p <- ggplot(PE10_PRvsNx_ggplot2_r, aes(x=DB,y=Dist.max)) + geom_violin() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PE10_PRxNx_Random2_violin.pdf")
p
dev.off()

# wo MAX
PE10_PRvsNXr_ggplot2woMax <- rbind(data.frame(PE10_PRvsNx_CCLEBPr,DB="GO_BP",Experiment="CCLE"),
  data.frame(PE10_PRvsNx_CCLEMFr,DB="GO_MF",Experiment="CCLE"),
 data.frame(PE10_PRvsNx_CCLECCr,DB="GO_CC",Experiment="CCLE"),
  data.frame(PE10_PRvsNx_GTEXBPr,DB="GO_BP",Experiment="GTEX"),
 data.frame(PE10_PRvsNx_GTEXMFr,DB="GO_MF",Experiment="GTEX"),
  data.frame(PE10_PRvsNx_GTEXCCr,DB="GO_CC",Experiment="GTEX"),
    data.frame(PE10_PRvsNx_TCGABPr,DB="GO_BP",Experiment="TCGA"),
  data.frame(PE10_PRvsNx_TCGAMFr,DB="GO_MF",Experiment="TCGA"),
   data.frame(PE10_PRvsNx_TCGACCr,DB="GO_CC",Experiment="TCGA"))
PE10_PRvsNXr_ggplot2woMax$DB <- factor(PE10_PRvsNXr_ggplot2woMax$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(PE10_PRvsNXr_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[1:3]) + scale_fill_manual(values=brewer.pal(4,"Set1")[1:3])
pdf("PE10_PRxNXr_Random2iteration1_histogram_woMax.pdf")
p
dev.off()

##########################
# 10% PE1 x CAFA
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_CCLEBPr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_CCLEBPr <- tmp

  } else {

    PE10_PRvsCAFA_CCLEBPr <- rbind(PE10_PRvsCAFA_CCLEBPr, tmp)

  }

}
colnames(PE10_PRvsCAFA_CCLEBPr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_CCLEBPr$Dist <- as.numeric(paste(PE10_PRvsCAFA_CCLEBPr$Dist))
PE10_PRvsCAFA_CCLEBPr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLEBPr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_CCLEBPr_max[PE10_PRvsCAFA_CCLEBPr_max==-Inf] <- NA

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_CCLEMFr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_CCLEMFr <- tmp

  } else {

    PE10_PRvsCAFA_CCLEMFr <- rbind(PE10_PRvsCAFA_CCLEMFr, tmp)

  }

}
colnames(PE10_PRvsCAFA_CCLEMFr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_CCLEMFr$Dist <- as.numeric(paste(PE10_PRvsCAFA_CCLEMFr$Dist))
PE10_PRvsCAFA_CCLEMFr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLEMFr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_CCLEMFr_max[PE10_PRvsCAFA_CCLEMFr_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_CCLECCr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_CCLECCr <- tmp

  } else {

    PE10_PRvsCAFA_CCLECCr <- rbind(PE10_PRvsCAFA_CCLECCr, tmp)

  }

}
colnames(PE10_PRvsCAFA_CCLECCr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_CCLECCr$Dist <- as.numeric(paste(PE10_PRvsCAFA_CCLECCr$Dist))
PE10_PRvsCAFA_CCLECCr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLECCr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_CCLECCr_max[PE10_PRvsCAFA_CCLECCr_max==-Inf] <- NA

PE10_PRvsCAFA_CCLEr <- rbind(PE10_PRvsCAFA_CCLEBPr_max, PE10_PRvsCAFA_CCLEMFr_max, PE10_PRvsCAFA_CCLECCr_max)
PE10_PRvsCAFA_CCLEr[,1] <- sapply(paste(PE10_PRvsCAFA_CCLEr[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_CCLEr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLEr, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFA_CCLEr_max)[2] <- "Dist.max"
PE10_PRvsCAFA_CCLEr_max[PE10_PRvsCAFA_CCLEr_max==-Inf] <- NA
 
# TCGA
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_TCGABPr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_TCGABPr <- tmp

  } else {

    PE10_PRvsCAFA_TCGABPr <- rbind(PE10_PRvsCAFA_TCGABPr, tmp)

  }

}
colnames(PE10_PRvsCAFA_TCGABPr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_TCGABPr$Dist <- as.numeric(paste(PE10_PRvsCAFA_TCGABPr$Dist))
PE10_PRvsCAFA_TCGABPr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGABPr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_TCGABPr_max[PE10_PRvsCAFA_TCGABPr_max==-Inf] <- NA

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_TCGAMFr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_TCGAMFr <- tmp

  } else {

    PE10_PRvsCAFA_TCGAMFr <- rbind(PE10_PRvsCAFA_TCGAMFr, tmp)

  }

}
colnames(PE10_PRvsCAFA_TCGAMFr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_TCGAMFr$Dist <- as.numeric(paste(PE10_PRvsCAFA_TCGAMFr$Dist))
PE10_PRvsCAFA_TCGAMFr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGAMFr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_TCGAMFr_max[PE10_PRvsCAFA_TCGAMFr_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_TCGACCr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_TCGACCr <- tmp

  } else {

    PE10_PRvsCAFA_TCGACCr <- rbind(PE10_PRvsCAFA_TCGACCr, tmp)

  }

}
colnames(PE10_PRvsCAFA_TCGACCr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_TCGACCr$Dist <- as.numeric(paste(PE10_PRvsCAFA_TCGACCr$Dist))
PE10_PRvsCAFA_TCGACCr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGACCr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_TCGACCr_max[PE10_PRvsCAFA_TCGACCr_max==-Inf] <- NA

PE10_PRvsCAFA_TCGAr <- rbind(PE10_PRvsCAFA_TCGABPr_max,PE10_PRvsCAFA_TCGAMFr_max,PE10_PRvsCAFA_TCGACCr_max)
PE10_PRvsCAFA_TCGAr[,1] <- sapply(paste(PE10_PRvsCAFA_TCGAr[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_TCGAr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGAr, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFA_TCGAr_max)[2] <- "Dist.max"
PE10_PRvsCAFA_TCGAr_max[PE10_PRvsCAFA_TCGAr_max==-Inf] <- NA

# GTEX
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_GTEXBPr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_GTEXBPr <- tmp

  } else {

    PE10_PRvsCAFA_GTEXBPr <- rbind(PE10_PRvsCAFA_GTEXBPr, tmp)

  }

}
colnames(PE10_PRvsCAFA_GTEXBPr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_GTEXBPr$Dist <- as.numeric(paste(PE10_PRvsCAFA_GTEXBPr$Dist))
PE10_PRvsCAFA_GTEXBPr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEXBPr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_GTEXBPr_max[PE10_PRvsCAFA_GTEXBPr_max==-Inf] <- NA

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_GTEXMFr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_GTEXMFr <- tmp

  } else {

    PE10_PRvsCAFA_GTEXMFr <- rbind(PE10_PRvsCAFA_GTEXMFr, tmp)

  }
}
colnames(PE10_PRvsCAFA_GTEXMFr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_GTEXMFr$Dist <- as.numeric(paste(PE10_PRvsCAFA_GTEXMFr$Dist))
PE10_PRvsCAFA_GTEXMFr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEXMFr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_GTEXMFr_max[PE10_PRvsCAFA_GTEXMFr_max==-Inf] <- NA

Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX/ALLEATORY_GOS_2/",full.names=T)
PE10_PRvsCAFA_GTEXCCr <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_GTEXCCr <- tmp

  } else {

    PE10_PRvsCAFA_GTEXCCr <- rbind(PE10_PRvsCAFA_GTEXCCr, tmp)

  }

}
colnames(PE10_PRvsCAFA_GTEXCCr) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_GTEXCCr$Dist <- as.numeric(paste(PE10_PRvsCAFA_GTEXCCr$Dist))
PE10_PRvsCAFA_GTEXCCr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEXCCr[,c(1,5)], FUN=max, na.rm=T)
PE10_PRvsCAFA_GTEXCCr_max[PE10_PRvsCAFA_GTEXCCr_max==-Inf] <- NA

PE10_PRvsCAFA_GTEXr <- rbind(PE10_PRvsCAFA_GTEXBPr_max,PE10_PRvsCAFA_GTEXMFr_max,PE10_PRvsCAFA_GTEXCCr_max)
PE10_PRvsCAFA_GTEXr[,1] <- sapply(paste(PE10_PRvsCAFA_GTEXr[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_GTEXr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEXr, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFA_GTEXr_max)[2] <- "Dist.max"
PE10_PRvsCAFA_GTEXr_max[PE10_PRvsCAFA_GTEXr_max==-Inf] <- NA

## ALL DBs and Datasets
PE10_PRvsCAFAr <- rbind(PE10_PRvsCAFA_CCLEr_max,PE10_PRvsCAFA_TCGAr_max,PE10_PRvsCAFA_GTEXr_max)
PE10_PRvsCAFAr_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFAr, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFAr_max)[2] <- "Dist.max"
PE10_PRvsCAFAr_max[PE10_PRvsCAFAr_max==-Inf] <- NA

PE10_PRvsCAFAr_ggplot2 <- rbind(data.frame(PE10_PRvsCAFA_CCLEBPr_max,DB="GO_BP",Experiment="CCLE"),
  data.frame(PE10_PRvsCAFA_CCLEMFr_max,DB="GO_MF",Experiment="CCLE"),
 data.frame(PE10_PRvsCAFA_CCLECCr_max,DB="GO_CC",Experiment="CCLE"),
 data.frame(PE10_PRvsCAFA_CCLEr_max,DB="All",Experiment="CCLE"),
 data.frame(PE10_PRvsCAFA_GTEXBPr_max,DB="GO_BP",Experiment="GTEX"),
 data.frame(PE10_PRvsCAFA_GTEXMFr_max,DB="GO_MF",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_GTEXCCr_max,DB="GO_CC",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_GTEXr_max,DB="All",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_TCGABPr_max,DB="GO_BP",Experiment="TCGA"),
  data.frame(PE10_PRvsCAFA_TCGAMFr_max,DB="GO_MF",Experiment="TCGA"),
   data.frame(PE10_PRvsCAFA_TCGACCr_max,DB="GO_CC",Experiment="TCGA"),
  data.frame(PE10_PRvsCAFA_TCGAr_max,DB="All",Experiment="TCGA"))
PE10_PRvsCAFAr_ggplot2$DB <- factor(PE10_PRvsCAFAr_ggplot2$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(PE10_PRvsCAFAr_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")) + scale_fill_manual(values=brewer.pal(4,"Set1"))
pdf("PE10_PRxCAFAr_Random2_histogram.pdf")
p
dev.off()
p <- ggplot(PE10_PRvsCAFAr_ggplot2, aes(x=DB,y=Dist.max)) + geom_boxplot() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PE10_PRxCAFAr_Random2_boxplot.pdf")
p
dev.off()
p <- ggplot(PE10_PRvsCAFAr_ggplot2, aes(x=DB,y=Dist.max)) + geom_violin() + ylab("max(Distance)") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(Experiment ~ .)
pdf("PE10_PRxCAFAr_Random2_violin.pdf")
p
dev.off()

# wo MAX
PE10_PRvsCAFAr_ggplot2woMax <- rbind(data.frame(PE10_PRvsCAFA_CCLEBPr,DB="GO_BP",Experiment="CCLE"),
  data.frame(PE10_PRvsCAFA_CCLEMFr,DB="GO_MF",Experiment="CCLE"),
 data.frame(PE10_PRvsCAFA_CCLECCr,DB="GO_CC",Experiment="CCLE"),
  data.frame(PE10_PRvsCAFA_GTEXBPr,DB="GO_BP",Experiment="GTEX"),
 data.frame(PE10_PRvsCAFA_GTEXMFr,DB="GO_MF",Experiment="GTEX"),
  data.frame(PE10_PRvsCAFA_GTEXCCr,DB="GO_CC",Experiment="GTEX"),
    data.frame(PE10_PRvsCAFA_TCGABPr,DB="GO_BP",Experiment="TCGA"),
  data.frame(PE10_PRvsCAFA_TCGAMFr,DB="GO_MF",Experiment="TCGA"),
   data.frame(PE10_PRvsCAFA_TCGACCr,DB="GO_CC",Experiment="TCGA"))
PE10_PRvsCAFAr_ggplot2woMax$DB <- factor(PE10_PRvsCAFAr_ggplot2woMax$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(PE10_PRvsCAFAr_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[1:3]) + scale_fill_manual(values=brewer.pal(4,"Set1")[1:3])
pdf("PE10_PRxCAFAr_Random2iteration1_histogram_woMax.pdf")
p
dev.off()


##############
# STATISTICS
p_cafa_NullDistrib_BP



stats <- data.frame(Total_uPE=c(nrow(PE10_PRvsNx_CCLEBP_max),nrow(PE10_PRvsNx_CCLEMF_max),nrow(PE10_PRvsNx_CCLECC_max),nrow(PE10_PRvsNx_CCLE_max), nrow(PE10_PRvsNx_TCGABP_max),nrow(PE10_PRvsNx_TCGAMF_max),nrow(PE10_PRvsNx_TCGACC_max),nrow(PE10_PRvsNx_TCGA_max), nrow(PE10_PRvsNx_GTEXBP_max),nrow(PE10_PRvsNx_GTEXMF_max),nrow(PE10_PRvsNx_GTEXCC_max),nrow(PE10_PRvsNx_GTEX_max), nrow(PE10_PRvsNx_max),
  nrow(PE10_PRvsCAFA_CCLEBP_max),nrow(PE10_PRvsCAFA_CCLEMF_max),nrow(PE10_PRvsCAFA_CCLECC_max),nrow(PE10_PRvsCAFA_CCLE_max), nrow(PE10_PRvsCAFA_TCGABP_max),nrow(PE10_PRvsCAFA_TCGAMF_max),nrow(PE10_PRvsCAFA_TCGACC_max),nrow(PE10_PRvsCAFA_TCGA_max), nrow(PE10_PRvsCAFA_GTEXBP_max),nrow(PE10_PRvsCAFA_GTEXMF_max),nrow(PE10_PRvsCAFA_GTEXCC_max),nrow(PE10_PRvsCAFA_GTEX_max), nrow(PE10_PRvsCAFA_max)),
  TP=c(sum(PE10_PRvsNx_CCLEBP_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_CCLEMF_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_CCLECC_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_CCLE_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_TCGABP_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_TCGAMF_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_TCGACC_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_TCGA_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_GTEXBP_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_GTEXMF_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_GTEXCC_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_GTEX_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_max[,2]==1,na.rm=T)
    ,sum(PE10_PRvsCAFA_CCLEBP_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_CCLEMF_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_CCLECC_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_CCLE_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_TCGABP_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_TCGAMF_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_TCGACC_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_TCGA_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_GTEXBP_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_GTEXMF_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_GTEXCC_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_GTEX_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_max[,2]==1,na.rm=T)), 
  DB=c("BP","MF","CC","All","BP","MF","CC","All","BP","MF","CC","All","All" ,"BP","MF","CC","All","BP","MF","CC","All","BP","MF","CC","All","All"), 
  Experiment=c("CCLE","CCLE","CCLE","CCLE","TCGA","TCGA","TCGA","TCGA","GTEX","GTEX","GTEX","GTEX","All" ,"CCLE","CCLE","CCLE","CCLE","TCGA","TCGA","TCGA","TCGA","GTEX","GTEX","GTEX","GTEX","All"), 
  REF=c("NX","NX","NX","NX","NX","NX","NX","NX","NX","NX","NX","NX","NX" ,"CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA"))

stats <- rbind(stats,data.frame(Total_uPE=c(nrow(PE10_PRvsNx_CCLEBPr_max),nrow(PE10_PRvsNx_CCLEMFr_max),nrow(PE10_PRvsNx_CCLECCr_max),nrow(PE10_PRvsNx_CCLEr_max), nrow(PE10_PRvsNx_TCGABPr_max),nrow(PE10_PRvsNx_TCGAMFr_max),nrow(PE10_PRvsNx_TCGACCr_max),nrow(PE10_PRvsNx_TCGAr_max), nrow(PE10_PRvsNx_GTEXBPr_max),nrow(PE10_PRvsNx_GTEXMFr_max),nrow(PE10_PRvsNx_GTEXCCr_max),nrow(PE10_PRvsNx_GTEXr_max), nrow(PE10_PRvsNxr_max),
  nrow(PE10_PRvsCAFA_CCLEBPr_max),nrow(PE10_PRvsCAFA_CCLEMFr_max),nrow(PE10_PRvsCAFA_CCLECCr_max),nrow(PE10_PRvsCAFA_CCLEr_max), nrow(PE10_PRvsCAFA_TCGABPr_max),nrow(PE10_PRvsCAFA_TCGAMFr_max),nrow(PE10_PRvsCAFA_TCGACCr_max),nrow(PE10_PRvsCAFA_TCGAr_max), nrow(PE10_PRvsCAFA_GTEXBPr_max),nrow(PE10_PRvsCAFA_GTEXMFr_max),nrow(PE10_PRvsCAFA_GTEXCCr_max),nrow(PE10_PRvsCAFA_GTEXr_max), nrow(PE10_PRvsCAFAr_max)),
  TP=c(sum(PE10_PRvsNx_CCLEBPr_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_CCLEMFr_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_CCLECCr_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_CCLEr_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_TCGABPr_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_TCGAMFr_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_TCGACCr_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_TCGAr_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_GTEXBPr_max[,2]==1,na.rm=T), sum(PE10_PRvsNx_GTEXMFr_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_GTEXCCr_max[,2]==1,na.rm=T),sum(PE10_PRvsNx_GTEXr_max[,2]==1,na.rm=T), sum(PE10_PRvsNxr_max[,2]==1,na.rm=T)
    ,sum(PE10_PRvsCAFA_CCLEBPr_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_CCLEMFr_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_CCLECCr_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_CCLEr_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_TCGABPr_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_TCGAMFr_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_TCGACCr_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_TCGAr_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_GTEXBPr_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_GTEXMFr_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_GTEXCCr_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_GTEXr_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFAr_max[,2]==1,na.rm=T)), 
  DB=c("BP","MF","CC","All","BP","MF","CC","All","BP","MF","CC","All","All" ,"BP","MF","CC","All","BP","MF","CC","All","BP","MF","CC","All","All"), 
  Experiment=c("CCLE_R","CCLE_R","CCLE_R","CCLE_R","TCGA_R","TCGA_R","TCGA_R","TCGA_R","GTEX_R","GTEX_R","GTEX_R","GTEX_R","All_R" ,"CCLE_R","CCLE_R","CCLE_R","CCLE_R","TCGA_R","TCGA_R","TCGA_R","TCGA_R","GTEX_R","GTEX_R","GTEX_R","GTEX_R","All_R"), 
  REF=c("NX","NX","NX","NX","NX","NX","NX","NX","NX","NX","NX","NX","NX" ,"CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA")))

stats <- stats[,c("REF","Experiment","DB","Total_uPE","TP")]
stats$Sens=stats$TP/stats$Total_uPE
write.table(stats, "stats_benchmarking.txt", row.names=F, sep="\t", quote=F)

# Difference between experiments
common_protPred <- c(intersect(paste(unique(PE10_PRvsCAFA_CCLEBP_max[,1])), intersect(paste(unique(PE10_PRvsCAFA_TCGABP_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXBP_max[,1])))), intersect(paste(unique(PE10_PRvsCAFA_CCLEMF_max[,1])), intersect(paste(unique(PE10_PRvsCAFA_TCGAMF_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXMF_max[,1])))), intersect(paste(unique(PE10_PRvsCAFA_CCLECC_max[,1])), intersect(paste(unique(PE10_PRvsCAFA_TCGACC_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXCC_max[,1])))))
names(common_protPred) <- c(rep("BP",length(intersect(paste(unique(PE10_PRvsCAFA_CCLEBP_max[,1])), intersect(paste(unique(PE10_PRvsCAFA_TCGABP_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXBP_max[,1])))))),rep("MF",length(intersect(paste(unique(PE10_PRvsCAFA_CCLEMF_max[,1])), intersect(paste(unique(PE10_PRvsCAFA_TCGAMF_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXMF_max[,1])))))), rep("CC", length(intersect(paste(unique(PE10_PRvsCAFA_CCLECC_max[,1])), intersect(paste(unique(PE10_PRvsCAFA_TCGACC_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXCC_max[,1])))))))
ccle_protPred <- c(setdiff(paste(unique(PE10_PRvsCAFA_CCLEBP_max[,1])), union(paste(unique(PE10_PRvsCAFA_TCGABP_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXBP_max[,1])))), setdiff(paste(unique(PE10_PRvsCAFA_CCLEMF_max[,1])), union(paste(unique(PE10_PRvsCAFA_TCGAMF_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXMF_max[,1])))), setdiff(paste(unique(PE10_PRvsCAFA_CCLECC_max[,1])), union(paste(unique(PE10_PRvsCAFA_TCGACC_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXCC_max[,1])))))
tcga_protPred <- c(setdiff(paste(unique(PE10_PRvsCAFA_TCGABP_max[,1])), union(paste(unique(PE10_PRvsCAFA_CCLEBP_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXBP_max[,1])))), setdiff(paste(unique(PE10_PRvsCAFA_TCGAMF_max[,1])), union(paste(unique(PE10_PRvsCAFA_CCLEMF_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXMF_max[,1])))), setdiff(paste(unique(PE10_PRvsCAFA_TCGACC_max[,1])), union(paste(unique(PE10_PRvsCAFA_CCLECC_max[,1])), paste(unique(PE10_PRvsCAFA_GTEXCC_max[,1])))))
gtex_protPred <- c(setdiff(paste(unique(PE10_PRvsCAFA_GTEXBP_max[,1])), union(paste(unique(PE10_PRvsCAFA_CCLEBP_max[,1])), paste(unique(PE10_PRvsCAFA_TCGABP_max[,1])))), setdiff(paste(unique(PE10_PRvsCAFA_GTEXMF_max[,1])), union(paste(unique(PE10_PRvsCAFA_CCLEMF_max[,1])), paste(unique(PE10_PRvsCAFA_TCGAMF_max[,1])))), setdiff(paste(unique(PE10_PRvsCAFA_GTEXCC_max[,1])), union(paste(unique(PE10_PRvsCAFA_CCLECC_max[,1])), paste(unique(PE10_PRvsCAFA_TCGACC_max[,1])))))
union(setdiff(paste(unique(PE10_PRvsCAFA_CCLECC_max[,1])),paste(common_protPred[names(common_protPred) %in% "CC"])), union(setdiff(paste(unique(PE10_PRvsCAFA_TCGACC_max[,1])),paste(common_protPred[names(common_protPred) %in% "CC"])), setdiff(paste(unique(PE10_PRvsCAFA_GTEXCC_max[,1])),paste(common_protPred[names(common_protPred) %in% "CC"]))) # 26 CC en 2 experimentos y 1401 en los 3

response_gg <- data.frame(Exp="CCLE",DB="BP", Distance=PE10_PRvsCAFA_CCLEBP[!is.na(PE10_PRvsCAFA_CCLEBP$Dist),"Dist"])
response_gg <- rbind(response_gg, data.frame(Exp="TCGA",DB="BP", Distance=PE10_PRvsCAFA_TCGABP[!is.na(PE10_PRvsCAFA_TCGABP$Dist),"Dist"]))
response_gg <- rbind(response_gg,data.frame(Exp="GTEX",DB="BP", Distance=PE10_PRvsCAFA_GTEXBP[!is.na(PE10_PRvsCAFA_GTEXBP$Dist),"Dist"]))
response_gg <- rbind(response_gg,data.frame(Exp="CCLE",DB="MF", Distance=PE10_PRvsCAFA_CCLEMF[!is.na(PE10_PRvsCAFA_CCLEMF$Dist),"Dist"]))
response_gg <- rbind(response_gg,data.frame(Exp="TCGA",DB="MF", Distance=PE10_PRvsCAFA_TCGAMF[!is.na(PE10_PRvsCAFA_TCGAMF$Dist),"Dist"]))
response_gg <- rbind(response_gg,data.frame(Exp="GTEX",DB="MF", Distance=PE10_PRvsCAFA_GTEXMF[!is.na(PE10_PRvsCAFA_GTEXMF$Dist),"Dist"]))
response_gg <- rbind(response_gg,data.frame(Exp="CCLE",DB="CC", Distance=PE10_PRvsCAFA_CCLECC[!is.na(PE10_PRvsCAFA_CCLECC$Dist),"Dist"]))
response_gg <- rbind(response_gg,data.frame(Exp="TCGA",DB="CC", Distance=PE10_PRvsCAFA_TCGACC[!is.na(PE10_PRvsCAFA_TCGACC$Dist),"Dist"]))
response_gg <- rbind(response_gg,data.frame(Exp="GTEX",DB="CC", Distance=PE10_PRvsCAFA_GTEXCC[!is.na(PE10_PRvsCAFA_GTEXCC$Dist),"Dist"]))
p <- ggplot(response_gg, aes(x=Exp,y=Distance)) + geom_boxplot() + ylab("Similarity") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(DB ~ .)
pdf("ExperimentEvaluation_boxplot.pdf")
p
dev.off()
p <- ggplot(response_gg, aes(x=Distance, col=Exp)) + geom_density() + facet_grid(DB ~ .) + theme_bw()
pdf("ExperimentEvaluation_density.pdf")
p
dev.off()

# Con el máximo
responseMax_gg <- data.frame(Exp="CCLE",DB="BP", Distance=PE10_PRvsCAFA_CCLEBP_max[!is.na(PE10_PRvsCAFA_CCLEBP_max$Dist.max),"Dist.max"])
responseMax_gg <- rbind(responseMax_gg, data.frame(Exp="TCGA",DB="BP", Distance=PE10_PRvsCAFA_TCGABP_max[!is.na(PE10_PRvsCAFA_TCGABP_max$Dist.max),"Dist.max"]))
responseMax_gg <- rbind(responseMax_gg,data.frame(Exp="GTEX",DB="BP", Distance=PE10_PRvsCAFA_GTEXBP_max[!is.na(PE10_PRvsCAFA_GTEXBP_max$Dist.max),"Dist.max"]))
responseMax_gg <- rbind(responseMax_gg,data.frame(Exp="CCLE",DB="MF", Distance=PE10_PRvsCAFA_CCLEMF_max[!is.na(PE10_PRvsCAFA_CCLEMF_max$Dist.max),"Dist.max"]))
responseMax_gg <- rbind(responseMax_gg,data.frame(Exp="TCGA",DB="MF", Distance=PE10_PRvsCAFA_TCGAMF_max[!is.na(PE10_PRvsCAFA_TCGAMF_max$Dist.max),"Dist.max"]))
responseMax_gg <- rbind(responseMax_gg,data.frame(Exp="GTEX",DB="MF", Distance=PE10_PRvsCAFA_GTEXMF_max[!is.na(PE10_PRvsCAFA_GTEXMF_max$Dist.max),"Dist.max"]))
responseMax_gg <- rbind(responseMax_gg,data.frame(Exp="CCLE",DB="CC", Distance=PE10_PRvsCAFA_CCLECC_max[!is.na(PE10_PRvsCAFA_CCLECC_max$Dist.max),"Dist.max"]))
responseMax_gg <- rbind(responseMax_gg,data.frame(Exp="TCGA",DB="CC", Distance=PE10_PRvsCAFA_TCGACC_max[!is.na(PE10_PRvsCAFA_TCGACC_max$Dist.max),"Dist.max"]))
responseMax_gg <- rbind(responseMax_gg,data.frame(Exp="GTEX",DB="CC", Distance=PE10_PRvsCAFA_GTEXCC_max[!is.na(PE10_PRvsCAFA_GTEXCC_max$Dist.max),"Dist.max"]))
p <- ggplot(responseMax_gg, aes(x=Exp,y=Distance)) + geom_boxplot() + ylab("Similarity") + theme_bw() + theme(axis.title.x=element_blank()) + facet_grid(DB ~ .)
pdf("ExperimentEvaluation_MaxDist_boxplot.pdf")
p
dev.off()
p <- ggplot(responseMax_gg, aes(x=Distance, col=Exp)) + geom_density() + facet_grid(DB ~ .) + theme_bw()
pdf("ExperimentEvaluation_MaxDist_density.pdf")
p
dev.off()

rownames(PE10_PRvsCAFA_GTEXBP_max) <- PE10_PRvsCAFA_GTEXBP_max[,1]
rownames(PE10_PRvsCAFA_GTEXMF_max) <- PE10_PRvsCAFA_GTEXMF_max[,1]
rownames(PE10_PRvsCAFA_GTEXCC_max) <- PE10_PRvsCAFA_GTEXCC_max[,1]
rownames(PE10_PRvsCAFA_TCGABP_max) <- PE10_PRvsCAFA_TCGABP_max[,1]
rownames(PE10_PRvsCAFA_TCGAMF_max) <- PE10_PRvsCAFA_TCGAMF_max[,1]
rownames(PE10_PRvsCAFA_TCGACC_max) <- PE10_PRvsCAFA_TCGACC_max[,1]
rownames(PE10_PRvsCAFA_CCLEBP_max) <- PE10_PRvsCAFA_CCLEBP_max[,1]
rownames(PE10_PRvsCAFA_CCLEMF_max) <- PE10_PRvsCAFA_CCLEMF_max[,1]
rownames(PE10_PRvsCAFA_CCLECC_max) <- PE10_PRvsCAFA_CCLECC_max[,1]

concordance_BP <- cbind(PE10_PRvsCAFA_GTEXBP_max, PE10_PRvsCAFA_TCGABP_max[rownames(PE10_PRvsCAFA_GTEXBP_max),2], PE10_PRvsCAFA_CCLEBP_max[rownames(PE10_PRvsCAFA_GTEXBP_max),2])
concordance_BP[concordance_BP==-Inf] <- NA
concordance_BP$n1 <- apply(concordance_BP, 1, FUN=function(x) sum(as.numeric(paste(x[-1]))==1, na.rm=T))
colnames(concordance_BP) <- c("Prot","GTEX","TCGA","CCLE","n1")
concordance_MF <- cbind(PE10_PRvsCAFA_GTEXMF_max, PE10_PRvsCAFA_TCGAMF_max[rownames(PE10_PRvsCAFA_GTEXMF_max),2], PE10_PRvsCAFA_CCLEMF_max[rownames(PE10_PRvsCAFA_GTEXMF_max),2])
concordance_MF[concordance_MF==-Inf] <- NA
concordance_MF$n1 <- apply(concordance_MF, 1, FUN=function(x) sum(as.numeric(paste(x[-1]))==1, na.rm=T))
colnames(concordance_MF) <- c("Prot","GTEX","TCGA","CCLE","n1")
concordance_CC <- cbind(PE10_PRvsCAFA_GTEXCC_max, PE10_PRvsCAFA_TCGACC_max[rownames(PE10_PRvsCAFA_GTEXCC_max),2], PE10_PRvsCAFA_CCLECC_max[rownames(PE10_PRvsCAFA_GTEXCC_max),2])
concordance_CC[concordance_CC==-Inf] <- NA
concordance_CC$n1 <- apply(concordance_CC, 1, FUN=function(x) sum(as.numeric(paste(x[-1]))==1, na.rm=T))
colnames(concordance_CC) <- c("Prot","GTEX","TCGA","CCLE","n1")
# Me coincide números individuales pero no conjuntos de anotaciones
length(unique(c(gsub("_GOBP_output.txt","",paste(concordance_BP[which(concordance_BP[,"TCGA"]==1),1])),gsub("_GOMF_output.txt","",paste(concordance_MF[which(concordance_MF[,"TCGA"]==1),1])),gsub("_GOCC_output.txt","",paste(concordance_CC[which(concordance_CC[,"TCGA"]==1),1])))))

concordance_df <- rbind(cbind(DB="BP",concordance_BP), cbind(DB="MF",concordance_MF),cbind(DB="CC",concordance_CC))
concordance_df[,2] <- gsub("_GOCC_output.txt","",gsub("_GOMF_output.txt","",gsub("_GOBP_output.txt","",concordance_df[,2])))
write.table(concordance_df, "PredictionsConcordance.txt", row.names=F, sep="\t", quote=F)

library(reshape2)
library(ggplot2)
library(stringr)

uniqueProteins <- paste(unique(concordance_df[,2]))
uniqueProteinsPerExp <- matrix(, nrow = 0, ncol = length(uniqueProteins))
#each column is a protein
colnames(uniqueProteinsPerExp) <- paste(uniqueProteins)
exps <- c("CCLE","TCGA","GTEX")
for(i in 1:length(exps))
{
  uniqueProteinsPerExp <- rbind(uniqueProteinsPerExp, (paste(uniqueProteins) %in% paste(unique(concordance_df[which(concordance_df[,exps[i]]==1),2])))*1)
}
rownames(uniqueProteinsPerExp) <- paste(exps)
#each row is a protein
uniqueProteinsPerExp_t <- t(uniqueProteinsPerExp)
matrix_uniqueProteinsPerExp <- uniqueProteinsPerExp %*% uniqueProteinsPerExp_t
matrix_uniqueProteinsPerExp.m <- melt(matrix_uniqueProteinsPerExp)
matrix_uniqueProteinsPerExp.mp <- matrix_uniqueProteinsPerExp.m
matrix_uniqueProteinsPerExp.mp$value <- round(matrix_uniqueProteinsPerExp.m$value/1560,digits=3)
p <- ggplot(matrix_uniqueProteinsPerExp.mp, aes(Var1, Var2)) + geom_tile(aes(fill=value), color="black") + geom_text(aes(fill = value, label = value)) + scale_fill_gradient(low = "white", high = "steelblue", trans = "log", na.value="white") + theme(legend.title=element_blank(), legend.text=element_blank()) + scale_x_discrete(name = "Experiments", labels = function(x) str_wrap(x, width = 10)) + scale_y_discrete(name = "Experiments") + theme_bw()
pdf("SensPerExp_AllOntoloties.pdf", width=5, height=4, colormodel="rgb")
p
dev.off()

# acertadas por las 3
length(unique(concordance_df[which(concordance_df$n1==3),2])) # 661 + 19 = 680
length(setdiff(unique(concordance_df[which(concordance_df$n1==2),2]), unique(concordance_df[which(concordance_df$n1==3),2]))) # 362
length(setdiff(unique(concordance_df[which(concordance_df$n1==2 & concordance_df$TCGA==1 & concordance_df$CCLE==1),2]), unique(concordance_df[which(concordance_df$n1==3),2]))) # 362 (GTEX-TCGA: 173-12=161, GTEX-CCLE: 64-9=55, TCGA-CCLE: 144-17=127)
length(setdiff(unique(concordance_df[which(concordance_df$n1==1),2]), unique(concordance_df[which(concordance_df$n1>1),2]))) # 271
length(setdiff(unique(concordance_df[which(concordance_df$n1==1 & concordance_df$CCLE==1),2]), unique(concordance_df[which(concordance_df$n1>1),2]))) # 271 (CCLE: 84, GTEX: 105, TCGA: 109)
length(unique(concordance_df[concordance_df$n1>0,2])) #1294

# No salen los números porque hay que hacer max (puede haber varias anotaciones acertadas y cada una con distintas DB para la prote)
concordance_df_max <- summaryBy(.~Prot, concordance_df[,-c(1,6)], FUN=max, na.rm=T)
concordance_df_max[concordance_df_max==-Inf] <- NA
concordance_df_max$n1 <- apply(concordance_df_max[,-1], 1, FUN=function(x) sum(as.numeric(paste(x))==1, na.rm=T))

length(unique(concordance_df_max[which(concordance_df_max$n1==3),1])) # 707
length(setdiff(unique(concordance_df_max[which(concordance_df_max$n1==2),1]), unique(concordance_df_max[which(concordance_df_max$n1==3),1]))) # 342

length(setdiff(paste(unique(concordance_df_max[which(concordance_df_max$n1==2 & concordance_df_max$GTEX.max==1 & concordance_df_max$CCLE.max==1),1])), paste(unique(concordance_df_max[which(concordance_df_max$n1==3),1])))) # 362 (GTEX-TCGA: 162, GTEX-CCLE: 54, TCGA-CCLE: 126)

length(setdiff(unique(concordance_df_max[which(concordance_df_max$n1==1),1]), unique(concordance_df_max[which(concordance_df_max$n1>1),1]))) # 245
length(setdiff(unique(concordance_df_max[which(concordance_df_max$n1==1 & concordance_df_max$CCLE.max==1),1]), unique(concordance_df_max[which(concordance_df_max$n1>1),1]))) # 271 (CCLE: 70, GTEX: 87, TCGA: 88)
length(unique(concordance_df_max[concordance_df_max$n1>0,2])) #1294

uniqueProteins <- paste(unique(concordance_df_max[,1]))
uniqueProteinsPerExp <- matrix(, nrow = 0, ncol = length(uniqueProteins))
#each column is a protein
colnames(uniqueProteinsPerExp) <- paste(uniqueProteins)
exps <- c("CCLE.max","TCGA.max","GTEX.max")
for(i in 1:length(exps))
{
  uniqueProteinsPerExp <- rbind(uniqueProteinsPerExp, (paste(uniqueProteins) %in% paste(unique(concordance_df_max[which(concordance_df_max[,exps[i]]==1),1])))*1)
}
rownames(uniqueProteinsPerExp) <- paste(exps)
#each row is a protein
uniqueProteinsPerExp_t <- t(uniqueProteinsPerExp)
matrix_uniqueProteinsPerExp <- uniqueProteinsPerExp %*% uniqueProteinsPerExp_t
matrix_uniqueProteinsPerExp.m <- melt(matrix_uniqueProteinsPerExp)
matrix_uniqueProteinsPerExp.mp <- matrix_uniqueProteinsPerExp.m
matrix_uniqueProteinsPerExp.mp$value <- round(matrix_uniqueProteinsPerExp.m$value/1560,digits=4)*100
p <- ggplot(matrix_uniqueProteinsPerExp.mp, aes(Var1, Var2)) + geom_tile(aes(fill=value), color="black") + geom_text(aes(fill = value, label = value)) + scale_fill_gradient(low = "white", high = "steelblue", trans = "log", na.value="white") + theme(legend.title=element_blank(), legend.text=element_blank()) + scale_x_discrete(name = "Experiments", labels = function(x) str_wrap(x, width = 10)) + scale_y_discrete(name = "Experiments") + theme_bw()
pdf("SensPerExpMax_AllOntoloties.pdf", width=5, height=4, colormodel="rgb")
p
dev.off()


# histogramas
p1 <- ggplot(concordance_df_max, aes(x=GTEX.max, y=TCGA.max)) + geom_point(alpha = 0.3) + theme_bw() + xlab("GTEX") + ylab("TCGA")
p2 <- ggplot(concordance_df_max, aes(x=GTEX.max, y=CCLE.max)) + geom_point(alpha = 0.3) + theme_bw() + xlab("GTEX") + ylab("CCLE")
p3 <- ggplot(concordance_df_max, aes(x=TCGA.max, y=CCLE.max)) + geom_point(alpha = 0.3) + theme_bw() + xlab("TCGA") + ylab("CCLE")
pdf("CCLExTCGA_scatter.pdf", width=3, height=3)
p3
dev.off()

g <- ggplot(concordance_df_max)

###########################################
# 10% PE1 RANDOM: Random null distribution
###########################################

#########
# CAFA

Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE",full.names=T)
PE10_PRvsCAFA_CCLEBPrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOBP_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsCAFA_CCLEBPrnull5[[i]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsCAFA_CCLEBPrnull5[[count]] <- paste(tmp[,6])
    }
  }
}
PE10_PRvsCAFA_CCLEBPrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_CCLEBPrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_CCLEBPrnull2), GO_PR=unlist(PE10_PRvsCAFA_CCLEBPrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_CCLEBPrnull4)), Random_n=unlist(PE10_PRvsCAFA_CCLEBPrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA",full.names=T)
PE10_PRvsCAFA_TCGABPrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOBP_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsCAFA_TCGABPrnull5[[1]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsCAFA_TCGABPrnull5[[count]] <-  paste(tmp[,6])
    }
  }
}
PE10_PRvsCAFA_TCGABPrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_TCGABPrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_TCGABPrnull2), GO_PR=unlist(PE10_PRvsCAFA_TCGABPrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_TCGABPrnull4)), Random_n=unlist(PE10_PRvsCAFA_TCGABPrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX",full.names=T)
PE10_PRvsCAFA_GTEXBPrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOBP_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsCAFA_GTEXBPrnull5[[i]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsCAFA_GTEXBPrnull5[[count]] <- paste(tmp[,6])
    }
  }
}
PE10_PRvsCAFA_GTEXBPrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_GTEXBPrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_GTEXBPrnull2), GO_PR=unlist(PE10_PRvsCAFA_GTEXBPrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_GTEXBPrnull4)), Random_n=unlist(PE10_PRvsCAFA_GTEXBPrnull5))

Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE",full.names=T)
PE10_PRvsCAFA_CCLEMFrnull1 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOMF_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsCAFA_CCLEMFrnull1[[i]] <- paste(tmp[,1])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsCAFA_CCLEMFrnull1[[count]] <- paste(tmp[,1])
    }
  }
}
PE10_PRvsCAFA_CCLEMFrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_CCLEMFrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_CCLEMFrnull2), GO_PR=unlist(PE10_PRvsCAFA_CCLEMFrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_CCLEMFrnull4)), Random_n=unlist(PE10_PRvsCAFA_CCLEMFrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA",full.names=T)
PE10_PRvsCAFA_TCGAMFrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOMF_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsCAFA_TCGAMFrnull5[[1]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsCAFA_TCGAMFrnull5[[count]] <-  paste(tmp[,6])
    }
  }
}
PE10_PRvsCAFA_TCGAMFrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_TCGAMFrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_TCGAMFrnull2), GO_PR=unlist(PE10_PRvsCAFA_TCGAMFrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_TCGAMFrnull4)), Random_n=unlist(PE10_PRvsCAFA_TCGAMFrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX",full.names=T)
PE10_PRvsCAFA_GTEXMFrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOMF_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsCAFA_GTEXMFrnull5[[i]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsCAFA_GTEXMFrnull5[[count]] <- paste(tmp[,6])
    }
  }
}
PE10_PRvsCAFA_GTEXMFrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_GTEXMFrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_GTEXMFrnull2), GO_PR=unlist(PE10_PRvsCAFA_GTEXMFrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_GTEXMFrnull4)), Random_n=unlist(PE10_PRvsCAFA_GTEXMFrnull5))

Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/CCLE",full.names=T)
PE10_PRvsCAFA_CCLECCrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOCC_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsCAFA_CCLECCrnull5[[i]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsCAFA_CCLECCrnull5[[count]] <- paste(tmp[,6])
    }
  }
}
PE10_PRvsCAFA_CCLECCrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_CCLECCrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_CCLECCrnull2), GO_PR=unlist(PE10_PRvsCAFA_CCLECCrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_CCLECCrnull4)), Random_n=unlist(PE10_PRvsCAFA_CCLECCrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/TCGA",full.names=T)
PE10_PRvsCAFA_TCGACCrnull1 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOCC_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))
    if (i == 1 & j==1) {
      PE10_PRvsCAFA_TCGACCrnull1[[i]] <- paste(tmp[,1])
      count <- 1
    } else {
      count <- count+1
      PE10_PRvsCAFA_TCGACCrnull1[[count]] <- paste(tmp[,1])
    }
  }
}
PE10_PRvsCAFA_TCGACCrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_TCGACCrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_TCGACCrnull2), GO_PR=unlist(PE10_PRvsCAFA_TCGACCrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_TCGACCrnull4)), Random_n=unlist(PE10_PRvsCAFA_TCGACCrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_VS_CAFA/GTEX",full.names=T)
PE10_PRvsCAFA_GTEXCCrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOCC_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))
    if (i == 1 & j==1) {
      PE10_PRvsCAFA_GTEXCCrnull5[[i]] <- paste(tmp[,6])
      count <- 1
    } else {
      count <- count+1
      PE10_PRvsCAFA_GTEXCCrnull5[[count]] <- paste(tmp[,6])
    }
  }
}
PE10_PRvsCAFA_GTEXCCrnull_df <- data.frame(Prot=unlist(PE10_PRvsCAFA_TCGACCrnull1),GO_CAFA=unlist(PE10_PRvsCAFA_TCGACCrnull2), GO_PR=unlist(PE10_PRvsCAFA_TCGACCrnull3), Distance=as.numeric(unlist(PE10_PRvsCAFA_TCGACCrnull4)), Random_n=unlist(PE10_PRvsCAFA_TCGACCrnull5))


#########
# NEXPROT

Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE",full.names=T)
PE10_PRvsNX_CCLEBPrnull1 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOBP_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsNX_CCLEBPrnull1[[i]] <- paste(tmp[,1])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsNX_CCLEBPrnull1[[count]] <- paste(tmp[,1])
    }
  }
}
PE10_PRvsNX_CCLEBPrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_CCLEBPrnull1),GO_CAFA=unlist(PE10_PRvsNX_CCLEBPrnull2), GO_PR=unlist(PE10_PRvsNX_CCLEBPrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_CCLEBPrnull4)), Random_n=unlist(PE10_PRvsNX_CCLEBPrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA",full.names=T)
PE10_PRvsNX_TCGABPrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOBP_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsNX_TCGABPrnull5[[1]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsNX_TCGABPrnull5[[count]] <-  paste(tmp[,6])
    }
  }
}
PE10_PRvsNX_TCGABPrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_TCGABPrnull1),GO_CAFA=unlist(PE10_PRvsNX_TCGABPrnull2), GO_PR=unlist(PE10_PRvsNX_TCGABPrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_TCGABPrnull4)), Random_n=unlist(PE10_PRvsNX_TCGABPrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX",full.names=T)
PE10_PRvsNX_GTEXBPrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOBP_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsNX_GTEXBPrnull5[[i]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsNX_GTEXBPrnull5[[count]] <- paste(tmp[,6])
    }
  }
}
PE10_PRvsNX_GTEXBPrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_GTEXBPrnull1),GO_CAFA=unlist(PE10_PRvsNX_GTEXBPrnull2), GO_PR=unlist(PE10_PRvsNX_GTEXBPrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_GTEXBPrnull4)), Random_n=unlist(PE10_PRvsNX_GTEXBPrnull5))


Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE",full.names=T)
PE10_PRvsNX_CCLEMFrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOMF_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsNX_CCLEMFrnull5[[i]] <- paste(tmp[,6])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsNX_CCLEMFrnull5[[count]] <- paste(tmp[,6])
    }

  }
}
PE10_PRvsNX_CCLEMFrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_CCLEMFrnull1),GO_CAFA=unlist(PE10_PRvsNX_CCLEMFrnull2), GO_PR=unlist(PE10_PRvsNX_CCLEMFrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_CCLEMFrnull4)), Random_n=unlist(PE10_PRvsNX_CCLEMFrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA",full.names=T)
PE10_PRvsNX_TCGAMFrnull1 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOMF_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsNX_TCGAMFrnull1[[1]] <- paste(tmp[,1])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsNX_TCGAMFrnull1[[count]] <-  paste(tmp[,1])
    }
  }
}
PE10_PRvsNX_TCGAMFrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_TCGAMFrnull1),GO_CAFA=unlist(PE10_PRvsNX_TCGAMFrnull2), GO_PR=unlist(PE10_PRvsNX_TCGAMFrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_TCGAMFrnull4)), Random_n=unlist(PE10_PRvsNX_TCGAMFrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX",full.names=T)
PE10_PRvsNX_GTEXMFrnull1 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOMF_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsNX_GTEXMFrnull1[[i]] <- paste(tmp[,1])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsNX_GTEXMFrnull1[[count]] <- paste(tmp[,1])
    }
  }
}
PE10_PRvsNX_GTEXMFrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_GTEXMFrnull1),GO_CAFA=unlist(PE10_PRvsNX_GTEXMFrnull2), GO_PR=unlist(PE10_PRvsNX_GTEXMFrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_GTEXMFrnull4)), Random_n=unlist(PE10_PRvsNX_GTEXMFrnull5))

Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/CCLE",full.names=T)
PE10_PRvsNX_CCLECCrnull1 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOCC_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))

    if (i == 1 & j==1) {

      PE10_PRvsNX_CCLECCrnull1[[i]] <- paste(tmp[,1])
      count <- 1

    } else {
      count <- count+1
      PE10_PRvsNX_CCLECCrnull1[[count]] <- paste(tmp[,1])
    }
  }
}
PE10_PRvsNX_CCLECCrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_CCLECCrnull1),GO_CAFA=unlist(PE10_PRvsNX_CCLECCrnull2), GO_PR=unlist(PE10_PRvsNX_CCLECCrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_CCLECCrnull4)), Random_n=unlist(PE10_PRvsNX_CCLECCrnull5))

Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/TCGA",full.names=T)
PE10_PRvsNX_TCGACCrnull5 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOCC_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))
    if (i == 1 & j==1) {
      PE10_PRvsNX_TCGACCrnull5[[i]] <- paste(tmp[,6])
      count <- 1
    } else {
      count <- count+1
      PE10_PRvsNX_TCGACCrnull5[[count]] <- paste(tmp[,6])
    }
  }
}
PE10_PRvsNX_TCGACCrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_TCGACCrnull1),GO_CAFA=unlist(PE10_PRvsNX_TCGACCrnull2), GO_PR=unlist(PE10_PRvsNX_TCGACCrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_TCGACCrnull4)), Random_n=unlist(PE10_PRvsNX_TCGACCrnull5))
Files <- list.files(pattern="MUESTREO", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/PE1_selection_10_percent/10_percent_PE1selected_PR_vs_NEXTPROT/GTEX",full.names=T)
PE10_PRvsNX_GTEXCCrnull1 <- list()
for (i in 1:length(Files)) {    
  File_tmp <- list.files(pattern="GOCC_output.txt$", path=paste(Files[i],"/",sep=""),full.names=T)
  cat(i, "\n")
  for(j in 1:length(File_tmp)){
    tmp <- read.csv(file = File_tmp[j], header = FALSE, col.names=c("V2","V3","V4","V5"), sep = " ", fill = TRUE)
    tmp <- cbind(unlist(strsplit(File_tmp[j], "//"))[2],tmp,substr(Files[i],111,nchar(Files[i])))
    if (i == 1 & j==1) {
      PE10_PRvsNX_GTEXCCrnull1[[i]] <- paste(tmp[,1])
      count <- 1
    } else {
      count <- count+1
      PE10_PRvsNX_GTEXCCrnull1[[count]] <- paste(tmp[,1])
    }
  }
}
PE10_PRvsNX_GTEXCCrnull_df <- data.frame(Prot=unlist(PE10_PRvsNX_GTEXCCrnull1),GO_CAFA=unlist(PE10_PRvsNX_GTEXCCrnull2), GO_PR=unlist(PE10_PRvsNX_GTEXCCrnull3), Distance=as.numeric(unlist(PE10_PRvsNX_GTEXCCrnull4)), Random_n=unlist(PE10_PRvsNX_GTEXCCrnull5))

randomNulldistribution_ggplot2 <- data.frame(Prot=c(unlist(PE10_PRvsCAFA_CCLEMFrnull1),unlist(PE10_PRvsCAFA_TCGAMFrnull1),unlist(PE10_PRvsCAFA_GTEXMFrnull1),unlist(PE10_PRvsCAFA_CCLECCrnull1),unlist(PE10_PRvsCAFA_TCGACCrnull1),unlist(PE10_PRvsCAFA_GTEXCCrnull1)),Distance=c(as.numeric(unlist(PE10_PRvsCAFA_CCLEMFrnull4)),as.numeric(unlist(PE10_PRvsCAFA_TCGAMFrnull4)),as.numeric(unlist(PE10_PRvsCAFA_GTEXMFrnull4)),as.numeric(unlist(PE10_PRvsCAFA_CCLECCrnull4)),as.numeric(unlist(PE10_PRvsCAFA_TCGACCrnull4)),as.numeric(unlist(PE10_PRvsCAFA_GTEXCCrnull4))), Random_n=c(unlist(PE10_PRvsCAFA_CCLEMFrnull5),unlist(PE10_PRvsCAFA_TCGAMFrnull5),unlist(PE10_PRvsCAFA_GTEXMFrnull5),unlist(PE10_PRvsCAFA_CCLECCrnull5),unlist(PE10_PRvsCAFA_TCGACCrnull5),unlist(PE10_PRvsCAFA_GTEXCCrnull5)), Experiment=c(rep("CCLE", length(unlist(PE10_PRvsCAFA_CCLEMFrnull1))), rep("TCGA", length(unlist(PE10_PRvsCAFA_TCGAMFrnull1))), rep("GTEX", length(unlist(PE10_PRvsCAFA_GTEXMFrnull1))),rep("CCLE", length(unlist(PE10_PRvsCAFA_CCLECCrnull1))), rep("TCGA", length(unlist(PE10_PRvsCAFA_TCGACCrnull1))), rep("GTEX", length(unlist(PE10_PRvsCAFA_GTEXCCrnull1)))), DB=c(rep("GO_MF",length(unlist(PE10_PRvsCAFA_CCLEMFrnull1))+length(unlist(PE10_PRvsCAFA_TCGAMFrnull1))+length(unlist(PE10_PRvsCAFA_GTEXMFrnull1))),rep("GO_CC",length(unlist(PE10_PRvsCAFA_CCLECCrnull1))+length(unlist(PE10_PRvsCAFA_TCGACCrnull1))+length(unlist(PE10_PRvsCAFA_GTEXCCrnull1)))))
randomNulldistribution_ggplot2$DB <- factor(paste(randomNulldistribution_ggplot2$DB),levels=c("GO_MF","GO_CC"))
# save(list=c("randomNulldistribution_ggplot2"), file="PRxCAFA_randomNullDistribution.rda")
# save(list=c("PE10_PRvsCAFA_CCLEBPrnull1","PE10_PRvsCAFA_CCLEBPrnull2","PE10_PRvsCAFA_CCLEBPrnull3","PE10_PRvsCAFA_CCLEBPrnull4","PE10_PRvsCAFA_CCLEBPrnull5"), file="PRxCAFA_BP_randomNullDistribution.rda")
# Me paso a sevastopol
randomNulldistributionFinal_ggplot2 <- data.frame(Prot=c(paste(PR_CAFA_CCLE_GOBP[,3]), paste(PR_CAFA_TCGA_GOBP[,3]), paste(PR_CAFA_GTEX_GOBP[,3]),paste(randomNulldistribution_ggplot2[,1])), Distance=c(PR_CAFA_CCLE_GOBP[,4], PR_CAFA_TCGA_GOBP[,4],PR_CAFA_GTEX_GOBP[,4],randomNulldistribution_ggplot2[,2]), Random_n=c(paste(PR_CAFA_CCLE_GOBP[,2]), paste(PR_CAFA_TCGA_GOBP[,2]), paste(PR_CAFA_GTEX_GOBP[,2]), paste(randomNulldistribution_ggplot2[,3])), Experiment=c(rep("CCLE",nrow(PR_CAFA_CCLE_GOBP)), rep("TCGA",nrow(PR_CAFA_TCGA_GOBP)), rep("GTEX",nrow(PR_CAFA_GTEX_GOBP)),paste(randomNulldistribution_ggplot2[,4])), DB=c(rep("GO_BP",nrow(PR_CAFA_CCLE_GOBP)), rep("GO_BP",nrow(PR_CAFA_TCGA_GOBP)), rep("GO_BP",nrow(PR_CAFA_GTEX_GOBP)), paste(randomNulldistribution_ggplot2[,5])))
# save(list=c("randomNulldistributionFinal_ggplot2"), file="PRxCAFA_randomNullDistribution_final.rda")
randomNulldistributionFinal_ggplot2$DB <- factor(paste(randomNulldistributionFinal_ggplot2$DB),levels=c("GO_BP", "GO_MF","GO_CC"))
randomNulldistributionFinal_ggplot2b <- randomNulldistributionFinal_ggplot2[-which(is.na(randomNulldistributionFinal_ggplot2[,2])),]
randomNulldistributionFinal_ggplot2b_filter <- randomNulldistributionFinal_ggplot2b[!(randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Random_n %in% levels(randomNulldistributionFinal_ggplot2b$Random_n)[51:100]),]
p <- ggplot(randomNulldistributionFinal_ggplot2b_filter, aes(x=Distance, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB, scales="free_y") + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[1:3]) + scale_fill_manual(values=brewer.pal(4,"Set1")[1:3])
scales="free_y"
pdf("PE10_PRxCAFAr_Random2_histogram.pdf")
p
dev.off()


PRxCAFA_stats_p05 <- data.frame(Exp=c(rep("CCLE",3),rep("TCGA",3),rep("GTEX",3),rep("ALL",3)), DB=rep(c("GO_BP","GO_MF","GO_CC"), 4),Total=c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC",])))
PRxCAFA_stats_p05$Distance_th05 <- c(sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[1]*0.05)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF","Distance"],decreasing=T)[round(PRxCAFA_stats_p05$Total[2]*0.05)],sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[3]*0.05)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[4]*0.05)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[5]*0.05)],sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[6]*0.05)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[7]*0.05)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[8]*0.05)],sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[9]*0.05)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[10]*0.05)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[11]*0.05)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[12]*0.05)]) 
PRxCAFA_stats_p05$greaterThan_th05 <-  c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[1],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[2],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[3],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[4],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[5],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[6],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[7],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[8],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[9],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[10],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[11],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th[12],]))
PRxCAFA_stats_p05$p_th05 <- PRxCAFA_stats_p05$greaterThan_th05/PRxCAFA_stats_p05$Total
PRxCAFA_stats_p05$Distance_th01 <- c(sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[1]*0.01)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF","Distance"],decreasing=T)[round(PRxCAFA_stats_p05$Total[2]*0.01)],sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[3]*0.01)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[4]*0.01)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[5]*0.01)],sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[6]*0.01)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[7]*0.01)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[8]*0.01)],sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[9]*0.01)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[10]*0.01)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[11]*0.01)], sort(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC","Distance"], decreasing=T)[round(PRxCAFA_stats_p05$Total[12]*0.01)]) 
PRxCAFA_stats_p05$greaterThan_th01 <-  c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[1],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[2],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[3],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[4],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[5],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[6],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[7],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[8],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[9],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[10],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[11],]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>PRxCAFA_stats_p05$Distance_th01[12],]))
PRxCAFA_stats_p05$p_th01 <- PRxCAFA_stats_p05$greaterThan_th01/PRxCAFA_stats_p05$Total
PRxCAFA_stats_p05$Distance_th1 <- rep(1,nrow(PRxCAFA_stats_p05)) 
PRxCAFA_stats_p05$equal_th1 <-  c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance==1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance==1,]))
PRxCAFA_stats_p05$p_th1 <- PRxCAFA_stats_p05$equal_th1/PRxCAFA_stats_p05$Total
PRxCAFA_stats <- data.frame(Exp=rep(c(rep("CCLE",3),rep("TCGA",3),rep("GTEX",3),rep("ALL",3)),5), DB=rep(c("GO_BP","GO_MF","GO_CC"), 20), Distance_th=c(rep(0.1,12),rep(0.3,12),rep(0.5,12),rep(0.7,12),rep(0.9,12)), Total=rep(c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF",]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC",])),5))
PRxCAFA_stats$Number_gt <- 0
PRxCAFA_stats$Number_gt[PRxCAFA_stats$Distance_th==0.1] <- c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.1,]))
PRxCAFA_stats$Number_gt[PRxCAFA_stats$Distance_th==0.3] <- c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.3,]))
PRxCAFA_stats$Number_gt[PRxCAFA_stats$Distance_th==0.5] <- c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.5,]))
PRxCAFA_stats$Number_gt[PRxCAFA_stats$Distance_th==0.7] <- c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.7,]))
PRxCAFA_stats$Number_gt[PRxCAFA_stats$Distance_th==0.9] <- c(nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="CCLE" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="TCGA" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$Exp=="GTEX" & randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_BP" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_MF" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]),nrow(randomNulldistributionFinal_ggplot2b[randomNulldistributionFinal_ggplot2b$DB=="GO_CC" & randomNulldistributionFinal_ggplot2b$Distance>0.9,]))

p_cafa_NullDistrib_BP <- c(PR_CAFA_CCLE_GOBP[,4], PR_CAFA_TCGA_GOBP[,4],PR_CAFA_GTEX_GOBP[,4])
p_cafa_NullDistrib_BP <- p_cafa_NullDistrib_BP[!is.na(p_cafa_NullDistrib_BP)]
p_cafa_NullDistrib_BP <- sort(p_cafa_NullDistrib_BP)
p_cafa_NullDistrib_BP_p <- sapply(c(0,unique(p_cafa_NullDistrib_BP)), FUN=function(x) sum(p_cafa_NullDistrib_BP>as.numeric(paste(x)))/length(p_cafa_NullDistrib_BP))
p_cafa_NullDistrib_BP_pDF <- data.frame(Distance=c(0,unique(p_cafa_NullDistrib_BP)), p=p_cafa_NullDistrib_BP_p)
rownames(p_cafa_NullDistrib_BP_pDF) <- p_cafa_NullDistrib_BP_pDF[,1]
p_cafa_NullDistrib_BP_pDF <- p_cafa_NullDistrib_BP_pDF[-1,]
# No están en la distribución nula los calculo a mano: "0.934" "0.935" "0.952" "0.905"
p_cafa_NullDistrib_BP_pDF <- rbind(p_cafa_NullDistrib_BP_pDF, "0.934"=c(0.934,(sum(p_cafa_NullDistrib_BP>0.934))/length(p_cafa_NullDistrib_BP)), "0.935"=c(0.935,(sum(p_cafa_NullDistrib_BP>0.935))/length(p_cafa_NullDistrib_BP)), "0.952"=c(0.952,(sum(p_cafa_NullDistrib_BP>0.952))/length(p_cafa_NullDistrib_BP)), "0.905"=c(0.905,(sum(p_cafa_NullDistrib_BP>0.905))/length(p_cafa_NullDistrib_BP)))

p_cafa_NullDistrib_MF <- randomNulldistribution_ggplot2$Distance[randomNulldistribution_ggplot2$DB=="GO_MF"]
p_cafa_NullDistrib_MF <- p_cafa_NullDistrib_MF[!is.na(p_cafa_NullDistrib_MF)]
p_cafa_NullDistrib_MF <- sort(p_cafa_NullDistrib_MF)
p_cafa_NullDistrib_MF_p <- sapply(c(0,unique(p_cafa_NullDistrib_MF)), FUN=function(x) sum(p_cafa_NullDistrib_MF>as.numeric(paste(x)))/length(p_cafa_NullDistrib_MF))
p_cafa_NullDistrib_MF_pDF <- data.frame(Distance=c(0,unique(p_cafa_NullDistrib_MF)), p=p_cafa_NullDistrib_MF_p)
rownames(p_cafa_NullDistrib_MF_pDF) <- p_cafa_NullDistrib_MF_pDF[,1]
p_cafa_NullDistrib_CC <- randomNulldistribution_ggplot2$Distance[randomNulldistribution_ggplot2$DB=="GO_CC"]
p_cafa_NullDistrib_CC <- p_cafa_NullDistrib_CC[!is.na(p_cafa_NullDistrib_CC)]
p_cafa_NullDistrib_CC <- sort(p_cafa_NullDistrib_CC)
p_cafa_NullDistrib_CC_p <- sapply(c(0,unique(p_cafa_NullDistrib_CC)), FUN=function(x) sum(p_cafa_NullDistrib_CC>as.numeric(paste(x)))/length(p_cafa_NullDistrib_CC))
p_cafa_NullDistrib_CC_pDF <- data.frame(Distance=c(0,unique(p_cafa_NullDistrib_CC)), p=p_cafa_NullDistrib_CC_p)
rownames(p_cafa_NullDistrib_CC_pDF) <- p_cafa_NullDistrib_CC_pDF[,1]

PE10_PRvsCAFA_CCLEBP$p <- p_cafa_NullDistrib_BP_pDF[paste(PE10_PRvsCAFA_CCLEBP$Dist), "p"]
PE10_PRvsCAFA_CCLEBP$fdr <- p.adjust(PE10_PRvsCAFA_CCLEBP$p, method="fdr")
PE10_PRvsCAFA_CCLEBP_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_CCLEBP[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_CCLEBP_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_CCLEBP[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_CCLEMF$p <- p_cafa_NullDistrib_MF_pDF[paste(PE10_PRvsCAFA_CCLEMF$Dist), "p"]
PE10_PRvsCAFA_CCLEMF$fdr <- p.adjust(PE10_PRvsCAFA_CCLEMF$p, method="fdr")
PE10_PRvsCAFA_CCLEMF_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_CCLEMF[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_CCLEMF_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_CCLEMF[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_CCLECC$p <- p_cafa_NullDistrib_CC_pDF[paste(PE10_PRvsCAFA_CCLECC$Dist), "p"]
PE10_PRvsCAFA_CCLECC$fdr <- p.adjust(PE10_PRvsCAFA_CCLECC$p, method="fdr")
PE10_PRvsCAFA_CCLECC_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_CCLECC[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_CCLECC_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_CCLECC[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_CCLE_max2 <- rbind(PE10_PRvsCAFA_CCLEBP_max2,PE10_PRvsCAFA_CCLEMF_max2,PE10_PRvsCAFA_CCLECC_max2)
PE10_PRvsCAFA_CCLE_max2[,1] <- sapply(paste(PE10_PRvsCAFA_CCLE_max2[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_CCLE_max2 <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLE_max2, FUN=min, na.rm=T)
colnames(PE10_PRvsCAFA_CCLE_max2)[2] <- "fdr.min"
PE10_PRvsCAFA_CCLE_max3 <- rbind(PE10_PRvsCAFA_CCLEBP_max3,PE10_PRvsCAFA_CCLEMF_max3,PE10_PRvsCAFA_CCLECC_max3)
PE10_PRvsCAFA_CCLE_max3[,1] <- sapply(paste(PE10_PRvsCAFA_CCLE_max3[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_CCLE_max3 <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_CCLE_max3, FUN=min, na.rm=T)
colnames(PE10_PRvsCAFA_CCLE_max3)[2] <- "p.min"
PE10_PRvsCAFA_TCGABP$p <- p_cafa_NullDistrib_BP_pDF[paste(PE10_PRvsCAFA_TCGABP$Dist), "p"]
PE10_PRvsCAFA_TCGABP$fdr <- p.adjust(PE10_PRvsCAFA_TCGABP$p, method="fdr")
PE10_PRvsCAFA_TCGABP_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_TCGABP[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_TCGABP_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_TCGABP[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_TCGAMF$p <- p_cafa_NullDistrib_MF_pDF[paste(PE10_PRvsCAFA_TCGAMF$Dist), "p"]
PE10_PRvsCAFA_TCGAMF$fdr <- p.adjust(PE10_PRvsCAFA_TCGAMF$p, method="fdr")
PE10_PRvsCAFA_TCGAMF_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_TCGAMF[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_TCGAMF_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_TCGAMF[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_TCGACC$p <- p_cafa_NullDistrib_CC_pDF[paste(PE10_PRvsCAFA_TCGACC$Dist), "p"]
PE10_PRvsCAFA_TCGACC$fdr <- p.adjust(PE10_PRvsCAFA_TCGACC$p, method="fdr")
PE10_PRvsCAFA_TCGACC_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_TCGACC[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_TCGACC_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_TCGACC[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_TCGA_max2 <- rbind(PE10_PRvsCAFA_TCGABP_max2,PE10_PRvsCAFA_TCGAMF_max2,PE10_PRvsCAFA_TCGACC_max2)
PE10_PRvsCAFA_TCGA_max2[,1] <- sapply(paste(PE10_PRvsCAFA_TCGA_max2[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_TCGA_max2 <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGA_max2, FUN=min, na.rm=T)
colnames(PE10_PRvsCAFA_TCGA_max2)[2] <- "fdr.min"
PE10_PRvsCAFA_TCGA_max3 <- rbind(PE10_PRvsCAFA_TCGABP_max3,PE10_PRvsCAFA_TCGAMF_max3,PE10_PRvsCAFA_TCGACC_max3)
PE10_PRvsCAFA_TCGA_max3[,1] <- sapply(paste(PE10_PRvsCAFA_TCGA_max3[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_TCGA_max3 <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_TCGA_max3, FUN=min, na.rm=T)
colnames(PE10_PRvsCAFA_TCGA_max3)[2] <- "p.min"

PE10_PRvsCAFA_GTEXBP$p <- p_cafa_NullDistrib_BP_pDF[paste(PE10_PRvsCAFA_GTEXBP$Dist), "p"]
PE10_PRvsCAFA_GTEXBP$fdr <- p.adjust(PE10_PRvsCAFA_GTEXBP$p, method="fdr")
PE10_PRvsCAFA_GTEXBP_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_GTEXBP[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_GTEXBP_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_GTEXBP[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_GTEXMF$p <- p_cafa_NullDistrib_MF_pDF[paste(PE10_PRvsCAFA_GTEXMF$Dist), "p"]
PE10_PRvsCAFA_GTEXMF$fdr <- p.adjust(PE10_PRvsCAFA_GTEXMF$p, method="fdr")
PE10_PRvsCAFA_GTEXMF_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_GTEXMF[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_GTEXMF_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_GTEXMF[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_GTEXCC$p <- p_cafa_NullDistrib_CC_pDF[paste(PE10_PRvsCAFA_GTEXCC$Dist), "p"]
PE10_PRvsCAFA_GTEXCC$fdr <- p.adjust(PE10_PRvsCAFA_GTEXCC$p, method="fdr")
PE10_PRvsCAFA_GTEXCC_max2 <- summaryBy(fdr ~ Prot,data=PE10_PRvsCAFA_GTEXCC[,c("Prot","fdr")], FUN=min, na.rm=T)
PE10_PRvsCAFA_GTEXCC_max3 <- summaryBy(p ~ Prot,data=PE10_PRvsCAFA_GTEXCC[,c("Prot","p")], FUN=min, na.rm=T)
PE10_PRvsCAFA_GTEX_max2 <- rbind(PE10_PRvsCAFA_GTEXBP_max2,PE10_PRvsCAFA_GTEXMF_max2,PE10_PRvsCAFA_GTEXCC_max2)
PE10_PRvsCAFA_GTEX_max2[,1] <- sapply(paste(PE10_PRvsCAFA_GTEX_max2[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_GTEX_max2 <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEX_max2, FUN=min, na.rm=T)
colnames(PE10_PRvsCAFA_GTEX_max2)[2] <- "fdr.min"
PE10_PRvsCAFA_GTEX_max3 <- rbind(PE10_PRvsCAFA_GTEXBP_max3,PE10_PRvsCAFA_GTEXMF_max3,PE10_PRvsCAFA_GTEXCC_max3)
PE10_PRvsCAFA_GTEX_max3[,1] <- sapply(paste(PE10_PRvsCAFA_GTEX_max3[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_GTEX_max3 <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_GTEX_max3, FUN=min, na.rm=T)
colnames(PE10_PRvsCAFA_GTEX_max3)[2] <- "p.min"
## ALL DBs and Datasets
PE10_PRvsCAFA_max2 <- rbind(PE10_PRvsCAFA_CCLE_max2,PE10_PRvsCAFA_TCGA_max2,PE10_PRvsCAFA_GTEX_max2)
PE10_PRvsCAFA_max2 <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_max2, FUN=min, na.rm=T)
colnames(PE10_PRvsCAFA_max2)[2] <- "fdr.min"
PE10_PRvsCAFA_max3 <- rbind(PE10_PRvsCAFA_CCLE_max3,PE10_PRvsCAFA_TCGA_max3,PE10_PRvsCAFA_GTEX_max3)
PE10_PRvsCAFA_max3 <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_max3, FUN=min, na.rm=T)
colnames(PE10_PRvsCAFA_max3)[2] <- "p.min"


statsFinal <- data.frame(Total_uPE=c(nrow(PE10_PRvsCAFA_CCLEBP_max),nrow(PE10_PRvsCAFA_CCLEMF_max),nrow(PE10_PRvsCAFA_CCLECC_max),nrow(PE10_PRvsCAFA_CCLE_max), nrow(PE10_PRvsCAFA_TCGABP_max),nrow(PE10_PRvsCAFA_TCGAMF_max),nrow(PE10_PRvsCAFA_TCGACC_max),nrow(PE10_PRvsCAFA_TCGA_max), nrow(PE10_PRvsCAFA_GTEXBP_max),nrow(PE10_PRvsCAFA_GTEXMF_max),nrow(PE10_PRvsCAFA_GTEXCC_max),nrow(PE10_PRvsCAFA_GTEX_max), nrow(PE10_PRvsCAFA_max)),

  TP_dist1=c(sum(PE10_PRvsCAFA_CCLEBP_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_CCLEMF_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_CCLECC_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_CCLE_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_TCGABP_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_TCGAMF_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_TCGACC_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_TCGA_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_GTEXBP_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_GTEXMF_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_GTEXCC_max[,2]==1,na.rm=T),sum(PE10_PRvsCAFA_GTEX_max[,2]==1,na.rm=T), sum(PE10_PRvsCAFA_max[,2]==1,na.rm=T)),

  TP_fdr05=c(sum(PE10_PRvsCAFA_CCLEBP_max2[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_CCLEMF_max2[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_CCLECC_max2[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_CCLE_max2[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_TCGABP_max2[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_TCGAMF_max2[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_TCGACC_max2[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_TCGA_max2[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_GTEXBP_max2[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_GTEXMF_max2[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_GTEXCC_max2[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_GTEX_max2[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_max2[,2]<0.05,na.rm=T)),

  TP_p05=c(sum(PE10_PRvsCAFA_CCLEBP_max3[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_CCLEMF_max3[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_CCLECC_max3[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_CCLE_max3[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_TCGABP_max3[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_TCGAMF_max3[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_TCGACC_max3[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_TCGA_max3[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_GTEXBP_max3[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_GTEXMF_max3[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_GTEXCC_max3[,2]<0.05,na.rm=T),sum(PE10_PRvsCAFA_GTEX_max3[,2]<0.05,na.rm=T), sum(PE10_PRvsCAFA_max3[,2]<0.05,na.rm=T)),

  DB=c("GO BP","GO MF","GO CC","All","GO BP","GO MF","GO CC","All","GO BP","GO MF","GO CC","All","All"), 
  Experiment=c("CCLE","CCLE","CCLE","CCLE","TCGA","TCGA","TCGA","TCGA","GTEX","GTEX","GTEX","GTEX","All"), 
  REF=c("CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA","CAFA"))
statsFinal <- statsFinal[,c(6,7,5,1:4)]
statsFinal$sensDist1 <- statsFinal[,5]/statsFinal[,4]
statsFinal$sensFDR05 <- statsFinal[,6]/statsFinal[,4]
statsFinal$sensP05 <- statsFinal[,7]/statsFinal[,4]
write.table(statsFinal, file="Final_stats_Benchmark.txt", sep="\t", quote=F)

statsFinal_ggplot2 <- rbind(data.frame(DB=statsFinal$DB, Proteins=statsFinal$Total_uPE, Number="Total", Experiment=statsFinal$Experiment), data.frame(DB=statsFinal$DB, Proteins=statsFinal$TP_dist1, Number="TP", Experiment=statsFinal$Experiment))
statsFinal_ggplot2$DB <- factor(statsFinal_ggplot2$DB, levels=c("GO BP", "GO MF", "GO CC", "All"))
statsFinal_ggplot2$Experiment <- factor(paste(statsFinal_ggplot2$Experiment), levels=c("CCLE","GTEX","TCGA","All"))
statsFinal_ggplot2$Number <- factor(statsFinal_ggplot2$Number, levels=c("Total","TP"))
p <- ggplot(statsFinal_ggplot2, aes(x=Experiment, y=Proteins, fill=Number)) + geom_bar(stat="identity", position=position_dodge()) + facet_grid(. ~ DB, scales="free") + theme_bw()
pdf("Benchmark_stats_dist1.pdf", 8, 8)
p
dev.off()
statsFinal_ggplot2 <- rbind(data.frame(DB=statsFinal$DB, Proteins=statsFinal$Total_uPE, Number="Total", Experiment=statsFinal$Experiment), data.frame(DB=statsFinal$DB, Proteins=statsFinal$TP_fdr05, Number="TP", Experiment=statsFinal$Experiment))
statsFinal_ggplot2$DB <- factor(statsFinal_ggplot2$DB, levels=c("GO BP", "GO MF", "GO CC", "All"))
statsFinal_ggplot2$Experiment <- factor(paste(statsFinal_ggplot2$Experiment), levels=c("CCLE","GTEX","TCGA","All"))
statsFinal_ggplot2$Number <- factor(statsFinal_ggplot2$Number, levels=c("Total","TP"))
p <- ggplot(statsFinal_ggplot2, aes(x=Experiment, y=Proteins, fill=Number)) + geom_bar(stat="identity", position=position_dodge()) + facet_grid(. ~ DB, scales="free") + theme_bw()
pdf("Benchmark_stats_fdr05.pdf", 8, 8)
p
dev.off()
statsFinal_ggplot2 <- rbind(data.frame(DB=statsFinal$DB, Proteins=statsFinal$Total_uPE, Number="Total", Experiment=statsFinal$Experiment), data.frame(DB=statsFinal$DB, Proteins=statsFinal$TP_p05, Number="TP", Experiment=statsFinal$Experiment))
statsFinal_ggplot2$DB <- factor(statsFinal_ggplot2$DB, levels=c("GO BP", "GO MF", "GO CC", "All"))
statsFinal_ggplot2$Experiment <- factor(paste(statsFinal_ggplot2$Experiment), levels=c("CCLE","GTEX","TCGA","All"))
statsFinal_ggplot2$Number <- factor(statsFinal_ggplot2$Number, levels=c("Total","TP"))
p <- ggplot(statsFinal_ggplot2, aes(x=Experiment, y=Proteins, fill=Number)) + geom_bar(stat="identity", position=position_dodge()) + facet_grid(. ~ DB, scales="free") + theme_bw()
pdf("Benchmark_stats_p05.pdf", 8, 8)
p
dev.off()


randomNulldistribution_nx_ggplot2 <- data.frame(Prot=c(unlist(PE10_PRvsNX_CCLEMFrnull1),unlist(PE10_PRvsNX_TCGAMFrnull1),unlist(PE10_PRvsNX_GTEXMFrnull1),unlist(PE10_PRvsNX_CCLECCrnull1),unlist(PE10_PRvsNX_TCGACCrnull1),unlist(PE10_PRvsNX_GTEXCCrnull1)),Distance=c(as.numeric(unlist(PE10_PRvsNX_CCLEMFrnull4)),as.numeric(unlist(PE10_PRvsNX_TCGAMFrnull4)),as.numeric(unlist(PE10_PRvsNX_GTEXMFrnull4)),as.numeric(unlist(PE10_PRvsNX_CCLECCrnull4)),as.numeric(unlist(PE10_PRvsNX_TCGACCrnull4)),as.numeric(unlist(PE10_PRvsNX_GTEXCCrnull4))), Random_n=c(unlist(PE10_PRvsNX_CCLEMFrnull5),unlist(PE10_PRvsNX_TCGAMFrnull5),unlist(PE10_PRvsNX_GTEXMFrnull5),unlist(PE10_PRvsNX_CCLECCrnull5),unlist(PE10_PRvsNX_TCGACCrnull5),unlist(PE10_PRvsNX_GTEXCCrnull5)), Experiment=c(rep("CCLE", length(unlist(PE10_PRvsNX_CCLEMFrnull1))), rep("TCGA", length(unlist(PE10_PRvsNX_TCGAMFrnull1))), rep("GTEX", length(unlist(PE10_PRvsNX_GTEXMFrnull1))),rep("CCLE", length(unlist(PE10_PRvsNX_CCLECCrnull1))), rep("TCGA", length(unlist(PE10_PRvsNX_TCGACCrnull1))), rep("GTEX", length(unlist(PE10_PRvsNX_GTEXCCrnull1)))), DB=c(rep("GO_MF",length(unlist(PE10_PRvsNX_CCLEMFrnull1))+length(unlist(PE10_PRvsNX_TCGAMFrnull1))+length(unlist(PE10_PRvsNX_GTEXMFrnull1))),rep("GO_CC",length(unlist(PE10_PRvsNX_CCLECCrnull1))+length(unlist(PE10_PRvsNX_TCGACCrnull1))+length(unlist(PE10_PRvsNX_GTEXCCrnull1)))))
randomNulldistribution_nx_ggplot2$DB <- factor(paste(randomNulldistribution_nx_ggplot2$DB),levels=c("GO_MF","GO_CC"))
# save(list=c("randomNulldistribution_nx_ggplot2"), file="PRxNX_randomNullDistribution.rda")
randomNulldistributionFinal_nx_ggplot2 <- data.frame(Prot=c(paste(PR_NX_CCLE_GOBP[,3]), paste(PR_NX_TCGA_GOBP[,3]), paste(PR_NX_GTEX_GOBP[,3]),paste(randomNulldistribution_nx_ggplot2[,1])), Distance=c(PR_NX_CCLE_GOBP[,4], PR_NX_TCGA_GOBP[,4],PR_NX_GTEX_GOBP[,4],randomNulldistribution_nx_ggplot2[,2]), Random_n=c(paste(PR_NX_CCLE_GOBP[,2]), paste(PR_NX_TCGA_GOBP[,2]), paste(PR_NX_GTEX_GOBP[,2]), paste(randomNulldistribution_nx_ggplot2[,3])), Experiment=c(rep("CCLE",nrow(PR_NX_CCLE_GOBP)), rep("TCGA",nrow(PR_NX_TCGA_GOBP)), rep("GTEX",nrow(PR_NX_GTEX_GOBP)),paste(randomNulldistribution_nx_ggplot2[,4])), DB=c(rep("GO_BP",nrow(PR_NX_CCLE_GOBP)), rep("GO_BP",nrow(PR_NX_TCGA_GOBP)), rep("GO_BP",nrow(PR_NX_GTEX_GOBP)), paste(randomNulldistribution_nx_ggplot2[,5])))

p <- ggplot(randomNulldistribution_nx_ggplot2, aes(x=Distance, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[2:3]) + scale_fill_manual(values=brewer.pal(4,"Set1")[2:3])
pdf("PE10_PRxNXrBP_Random2_histogram.pdf")
p
dev.off()


##########################
# PROTEÓMICA NCI60
##########################

##########################
# newPE1 x nextProt
# CCLE
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/NCI60_spectral_counts/NCI60_NEXTPROT/",full.names=T)
PE10_PRvsNx_NCI60BP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_NCI60BP <- tmp

  } else {

    PE10_PRvsNx_NCI60BP <- rbind(PE10_PRvsNx_NCI60BP, tmp)
  }
}
colnames(PE10_PRvsNx_NCI60BP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_NCI60BP$Dist <- as.numeric(paste(PE10_PRvsNx_NCI60BP$Dist))
PE10_PRvsNx_NCI60BP_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_NCI60BP[,c(1,5)], FUN=max, na.rm=T)

Files <- list.files(pattern="GOMF_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/NCI60_spectral_counts/NCI60_NEXTPROT/",full.names=T)
PE10_PRvsNx_NCI60MF <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_NCI60MF <- tmp

  } else {

    PE10_PRvsNx_NCI60MF <- rbind(PE10_PRvsNx_NCI60MF, tmp)

  }

}
colnames(PE10_PRvsNx_NCI60MF) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_NCI60MF$Dist <- as.numeric(paste(PE10_PRvsNx_NCI60MF$Dist))
PE10_PRvsNx_NCI60MF_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_NCI60MF[,c(1,5)], FUN=max, na.rm=T)
Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/NCI60_spectral_counts/NCI60_NEXTPROT/",full.names=T)
PE10_PRvsNx_NCI60CC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsNx_NCI60CC <- tmp

  } else {

    PE10_PRvsNx_NCI60CC <- rbind(PE10_PRvsNx_NCI60CC, tmp)

  }

}
colnames(PE10_PRvsNx_NCI60CC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsNx_NCI60CC$Dist <- as.numeric(paste(PE10_PRvsNx_NCI60CC$Dist))
PE10_PRvsNx_NCI60CC_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_NCI60CC[,c(1,5)], FUN=max, na.rm=T)

PE10_PRvsNx_NCI60 <- rbind(PE10_PRvsNx_NCI60BP_max,PE10_PRvsNx_NCI60MF_max,PE10_PRvsNx_NCI60CC_max)
PE10_PRvsNx_NCI60[,1] <- sapply(paste(PE10_PRvsNx_NCI60[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsNx_NCI60_max <- summaryBy(. ~ Prot, data=PE10_PRvsNx_NCI60, FUN=max, na.rm=T)

PE10_PRvsNx_NCI60_ggplot2 <- rbind(data.frame(PE10_PRvsNx_NCI60BP_max,DB="GO_BP",Experiment="NCI60"),
  data.frame(PE10_PRvsNx_NCI60MF_max,DB="GO_MF",Experiment="NCI60"),
 data.frame(PE10_PRvsNx_NCI60CC_max,DB="GO_CC",Experiment="NCI60"),
 data.frame(PE10_PRvsNx_NCI60_max,DB="All",Experiment="NCI60"))
PE10_PRvsNx_NCI60_ggplot2$DB <- factor(PE10_PRvsNx_NCI60_ggplot2$DB, levels=c("GO_BP","GO_MF","GO_CC","All"))

p <- ggplot(PE10_PRvsNx_NCI60_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() +
    scale_color_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1")
pdf("newPE_PRxNx_NCI60_histogram.pdf")
p
dev.off()

# wo MAX
PE10_PRvsNX_NCI60_ggplot2woMax <- rbind(data.frame(PE10_PRvsNx_NCI60BP,DB="GO_BP",Experiment="NCI60"),
  data.frame(PE10_PRvsNx_NCI60MF,DB="GO_MF",Experiment="NCI60"),
 data.frame(PE10_PRvsNx_NCI60CC,DB="GO_CC",Experiment="NCI60"))
PE10_PRvsNX_NCI60_ggplot2woMax$DB <- factor(PE10_PRvsNX_NCI60_ggplot2woMax$DB, levels=c("GO_BP","GO_MF","GO_CC"))

p <- ggplot(PE10_PRvsNX_NCI60_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[1:3]) + scale_fill_manual(values=brewer.pal(4,"Set1")[1:3])
pdf("PE10_PRxNX_NCI60_histogram_woMax.pdf")
p
dev.off()


##########################
# new PE1 x CAFA
# NCI60
Files <- list.files(pattern="GOBP_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/NCI60_spectral_counts/NCI60_CAFA/",full.names=T)
PE10_PRvsCAFA_NCI60BP <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_NCI60BP <- tmp

  } else {

    PE10_PRvsCAFA_NCI60BP <- rbind(PE10_PRvsCAFA_NCI60BP, tmp)

  }

}
colnames(PE10_PRvsCAFA_NCI60BP) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_NCI60BP$Dist <- as.numeric(paste(PE10_PRvsCAFA_NCI60BP$Dist))
PE10_PRvsCAFA_NCI60BP_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_NCI60BP[,c(1,5)], FUN=max, na.rm=T)
Files <- list.files(pattern="GOCC_output.txt$", path="../../../pepe/43_UPEFINDER_AGO20/Benchmarking/NCI60_spectral_counts/NCI60_CAFA/",full.names=T)
PE10_PRvsCAFA_NCI60CC <- data.frame()
for (i in 1:length(Files)) {    

  cat(i, "\n")

  tmp <- read.csv(file = Files[i], header = FALSE, sep = " ", col.names=c("V2","V3","V4","V5"), fill = TRUE)
  tmp <- cbind(unlist(strsplit(Files[i], "//"))[2],tmp)

  if (i == 1) {

    PE10_PRvsCAFA_NCI60CC <- tmp

  } else {

    PE10_PRvsCAFA_NCI60CC <- rbind(PE10_PRvsCAFA_NCI60CC, tmp)

  }

}
colnames(PE10_PRvsCAFA_NCI60CC) <- c("Prot","GO1","GO2","Ont","Dist")
PE10_PRvsCAFA_NCI60CC$Dist <- as.numeric(paste(PE10_PRvsCAFA_NCI60CC$Dist))
PE10_PRvsCAFA_NCI60CC_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_NCI60CC[,c(1,5)], FUN=max, na.rm=T)

PE10_PRvsCAFA_NCI60 <- rbind(PE10_PRvsCAFA_NCI60BP_max, PE10_PRvsCAFA_NCI60CC_max)
PE10_PRvsCAFA_NCI60[,1] <- sapply(paste(PE10_PRvsCAFA_NCI60[,1]), FUN=function(x) unlist(strsplit(x, "_GO"))[1])
PE10_PRvsCAFA_NCI60_max <- summaryBy(. ~ Prot, data=PE10_PRvsCAFA_NCI60, FUN=max, na.rm=T)
colnames(PE10_PRvsCAFA_NCI60_max)[2] <- "Dist.max"


PE10_PRvsCAFA_NCI60_ggplot2 <- rbind(data.frame(PE10_PRvsCAFA_NCI60BP_max,DB="GO_BP",Experiment="NCI60"),
  data.frame(PE10_PRvsCAFA_NCI60CC_max,DB="GO_CC",Experiment="NCI60"),
 data.frame(PE10_PRvsCAFA_NCI60_max,DB="All",Experiment="NCI60"))
PE10_PRvsCAFA_NCI60_ggplot2$DB <- factor(PE10_PRvsCAFA_NCI60_ggplot2$DB, levels=c("GO_BP","GO_CC","All"))
p <- ggplot(PE10_PRvsCAFA_NCI60_ggplot2, aes(x=Dist.max, col=DB, fill=DB)) + geom_histogram() + facet_grid(. ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[c(1,3,4)]) + scale_fill_manual(values=brewer.pal(4,"Set1")[c(1,3,4)])
pdf("newPE_PRxCAFA_NCI60_histogram.pdf")
p
dev.off()

# wo MAX
PE10_PRvsCAFA_NCI60_ggplot2woMax <- rbind(data.frame(PE10_PRvsCAFA_NCI60BP,DB="GO_BP",Experiment="NCI60"),
  data.frame(PE10_PRvsCAFA_NCI60CC,DB="GO_CC",Experiment="NCI60"))
PE10_PRvsCAFA_NCI60_ggplot2woMax$DB <- factor(PE10_PRvsCAFA_NCI60_ggplot2woMax$DB, levels=c("GO_BP","GO_CC"))

p <- ggplot(PE10_PRvsCAFA_NCI60_ggplot2woMax, aes(x=Dist, col=DB, fill=DB)) + geom_histogram() + facet_grid(Experiment ~ DB) + theme_bw() + scale_colour_manual(values=brewer.pal(4,"Set1")[c(1,3)]) + scale_fill_manual(values=brewer.pal(4,"Set1")[c(1,3)])
pdf("newPE_PRxCAFA_histogram_woMax.pdf")
p
dev.off()

