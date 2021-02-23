#!/usr/bin/Rscript

args=(commandArgs(TRUE))
projectDir <- args[1]
#databasename <- args[2]

library(Biostrings)
library(doBy)
library(reshape2)
library(ggplot2)
library(stringr)
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesVikv2.R")

mascot_rda_folder <- paste0(projectDir, "/Dat_Files/")


mascot_rdaFiles <- list.files(path = mascot_rda_folder, pattern = "*dataPSMMat_filterPSM.txt")

mascot_rdaFiles <- paste0(mascot_rda_folder,mascot_rdaFiles)

resultDir <- paste0(projectDir, "/Results/")
y <- paste0("mkdir -p ", resultDir)
system(y)
y <- NULL

for (i in 1:length(mascot_rdaFiles)) {
	data<-read.csv2(paste(mascot_rdaFiles[i]),header=T, sep="\t")
	cat(paste(mascot_rdaFiles[i]),"\n")
	if(i == 1) {
		results_mascot_TyD <- data
	} else {
		results_mascot_TyD <- rbind(results_mascot_TyD, data)
	}
}
save(results_mascot_TyD, file = paste0(resultDir, "results_mascot_TyD.rda"))
write.table(results_mascot_TyD,file = paste0(resultDir, "results_mascot_TyD.txt"), row.names=F, col.names=T, sep="\t", quote=F)

cat("Results Mascot:")
cat("\n")
cat("Range PSMFDR:", range(as.numeric(paste(results_mascot_TyD$psmFDR))))
cat("\n")
cat("Range pepFDR:", range(as.numeric(paste(results_mascot_TyD$pepFDR))))
cat("\n")
cat("Range protFDR:", range(as.numeric(paste(results_mascot_TyD$protFDR))))
cat("\n")
cat("Number of samples:", length(table(results_mascot_TyD$sample)))
cat("\n")
cat("Samples:", length(unique(results_mascot_TyD$sample)))
cat("\n")
cat("Number of target and decoy values:", table(results_mascot_TyD$database))
cat("\n")
cat("Example of target protein names:", paste(results_mascot_TyD[results_mascot_TyD$database == "T", "ProteinAccession"][1:5]))
cat("\n")
cat("Example of decoy protein names:", paste(results_mascot_TyD[results_mascot_TyD$database == "D", "ProteinAccession"][1:5]))
cat("\n")

save(results_mascot_TyD, file = paste0(resultDir, "results_mascot_TyD.rda"))

# remove decoy identifications

results_mascot <- results_mascot_TyD[results_mascot_TyD$database == "T",]

#results_mascot$PSM <- paste(results_mascot$datfile, results_mascot$Query, sep = "-")

save(results_mascot, file = paste0(resultDir, "results_mascot_onlyT.rda"))

for_crux<-results_mascot_TyD[,c("ScanString","score","psmFDR","PeptideSeq","ProteinAccession")]

names(for_crux)<-c("psm_id","score",	"q-value",	"sequence",	"proteinIds")

for_crux$psm_id<-str_replace_all(paste(for_crux$psm_id)," ","_")

results_sum<-aggregate(sequence ~ psm_id + proteinIds, data = unique(for_crux[,c("psm_id","sequence","proteinIds")]), FUN = paste, collapse = ",")
results_sum$count<-str_count(paste(results_sum$sequence),",")+1

tmp_res<-results_sum[,c("psm_id","count")]
tmp_res<-unique(tmp_res)

total_res<-merge(for_crux,tmp_res, by="psm_id", all=T)

names(total_res)<-c("psm_id", "score","q-value","sequence","protein id","distinct matches/spectrum")

write.table(total_res, file=paste0(resultDir, "results_mascot_TyD_TO_spectral_counting.txt"), row.names=F,col.names=T, sep="\t", quote=F)
###need(PSMId	score	q-value	posterior_error_prob	peptide	proteinIds)




cat("end of script","\n")

### out of R
#crux spectral-counts --verbosity 40 --threshold-type none --fileroot results_spec_counting_NSAF --output-dir /home/margaret/data/pepe/14_NCI60_uPI_FEB18/SPEC-COUNT-results/ --measure NSAF --protein-database /home/margaret/data/pepe/for_spec.fasta --spectrum-parser mstoolkit /home/margaret/data/pepe/14_NCI60_uPI_FEB18/Results/results_mascot_TyD_TO_spectral_counting.txt --overwrite T
