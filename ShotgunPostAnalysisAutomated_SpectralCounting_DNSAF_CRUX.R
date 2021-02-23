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


mascot_rdaFiles <- list.files(path = mascot_rda_folder, pattern = "*PSMFDR_filter.rda")

mascot_rdaFiles <- paste0(mascot_rda_folder,mascot_rdaFiles)

resultDir <- paste0(projectDir, "/Files_to_crux/")
y <- paste0("mkdir -p ", resultDir)
system(y)
y <- NULL

for (i in 1:length(mascot_rdaFiles)) {
	data<-get(load(mascot_rdaFiles[i]))
	cat(paste(mascot_rdaFiles[i]),"\n")



	cat("Results Mascot:")
	cat("\n")
	cat("Range PSMFDR:", range(as.numeric(paste(data$psmFDR))))
	cat("\n")
	cat("Range pepFDR:", range(as.numeric(paste(data$pepFDR))))
	cat("\n")
	cat("Range protFDR:", range(as.numeric(paste(data$protFDR))))
	cat("\n")
	cat("Number of samples:", length(table(data$sample)))
	cat("\n")
	cat("Samples:", length(unique(data$sample)))
	cat("\n")
	cat("Number of target and decoy values:", table(data$database))
	cat("\n")
	cat("Example of target protein names:", paste(data[data$database == "T", "ProteinAccession"][1:5]))
	cat("\n")
	cat("Example of decoy protein names:", paste(data[data$database == "D", "ProteinAccession"][1:5]))
	cat("\n")


# remove decoy identifications

#results_mascot <- results_mascot_TyD[results_mascot_TyD$database == "T",]

#results_mascot$PSM <- paste(results_mascot$datfile, results_mascot$Query, sep = "-")

#save(results_mascot, file = paste0(resultDir, "results_mascot_onlyT.rda"))

	for_crux<-data[,c("ScanString","score","psmFDR","PeptideSeq","ProteinAccession")]

	names(for_crux)<-c("psm_id","score",	"q-value",	"sequence",	"proteinIds")

	for_crux$psm_id<-str_replace_all(paste(for_crux$psm_id)," ","_")

	results_sum<-aggregate(sequence ~ psm_id + proteinIds, data = unique(for_crux[,c("psm_id","sequence","proteinIds")]), FUN = paste, collapse = ",")
	results_sum$count<-str_count(paste(results_sum$sequence),",")+1

	tmp_res<-results_sum[,c("psm_id","count")]
	tmp_res<-unique(tmp_res)

	total_res<-merge(for_crux,tmp_res, by="psm_id", all=T)

	names(total_res)<-c("psm_id", "score","q-value","sequence","protein id","distinct matches/spectrum")

	file_name<-basename(mascot_rdaFiles[i])
	name_file<-strsplit(file_name,"\\.")[[1]][1]

	new_name<-paste0(name_file,"_TO_crux.txt")
	new_exit<-paste0(resultDir,new_name)

	write.table(total_res, file=new_exit, row.names=F,col.names=T, sep="\t", quote=F)
###need(PSMId	score	q-value	posterior_error_prob	peptide	proteinIds)


}

	cat("end of script","\n")

### out of R
#crux spectral-counts --verbosity 40 --threshold-type none --fileroot results_spec_counting_NSAF --output-dir /home/margaret/data/pepe/14_NCI60_uPI_FEB18/SPEC-COUNT-results/ --measure NSAF --protein-database /home/margaret/data/pepe/for_spec.fasta --spectrum-parser mstoolkit /home/margaret/data/pepe/14_NCI60_uPI_FEB18/Results/results_mascot_TyD_TO_spectral_counting.txt --overwrite T
