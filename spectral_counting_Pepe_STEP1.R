
args=(commandArgs(TRUE))
experiment_Path<-args[1]   #### /home/margaret/data/pepe/14_NCI60_uPI_FEB18/
database<-args[2]          #### complete path to FASTA FILE (TARGET & DECOY IN ONE FILE)
method<-args[3]
log<-args[4]    #TRUE OR FALSE

library(doBy)
library(Biostrings)
library(stringr)

path_files<-paste0(experiment_Path,"Dat_Files/")
name_sample_files<-list.files(path=path_files, pattern="PSMFDR_filter.rda")
file_path<-paste0(path_files,name_sample_files)

fasta<-readAAStringSet(database)
lengths<-paste(width(fasta))
proteins<-paste(names(fasta))
db_target_nAA<-data.frame("lengths"=lengths, "ProteinACC"=proteins)
db_target_nAA$ProteinAccession<-lapply(strsplit(paste(db_target_nAA$ProteinACC)," "), "[", 1)
db_target_nAA$ProteinACC<-NULL

spectralCounting_Pepe_NSAF <- function(results_sample, db_target_nAA, log) {

	sample = paste(unique(results_sample$sample))
	cat("Calculating spectral counting ... ", sample," ... \n")
	tmp <- unique(results_sample[, c("datfile", "Query", "ProteinAccession")])
	tmp$count <- paste(tmp$datfile, tmp$Query, sep = "_")
	tmp <- tmp[, c("count", "ProteinAccession")]
	results_sample_Counts <- summaryBy(count ~ ProteinAccession, data = tmp, FUN = function(x) length(x))
	colnames(results_sample_Counts) <- c("ProteinAccession", "SpC")
	object_int<-merge(results_sample_Counts,db_target_nAA, by.x="ProteinAccession", by.y="ProteinAccession")

	#results_sample_Counts$length_protein <- db_target_nAA[paste(results_sample_Counts[,1])]
	results_sample_Counts<-object_int
	names(results_sample_Counts)<-c("ProteinAccession","SpC","length_protein")
	results_sample_Counts$NSAF <- (as.numeric(results_sample_Counts$SpC)/as.numeric(results_sample_Counts$length_protein))/sum(as.numeric(results_sample_Counts$SpC)/as.numeric(results_sample_Counts$length_protein))

	if(log=="T"){
		results_sample_Counts$NSAF<-log10(results_sample_Counts$NSAF)
	}

	results_sample_Counts$sample <- sample
	return(results_sample_Counts)
}

spectralCounting_Pepe_uNSAF<- function(results_sample, db_target_nAA, log){
	sample = paste(unique(results_sample$sample))
	cat("Calculating spectral counting ... ", sample," ... \n")
	tmp <- unique(results_sample[, c("datfile", "Query","PeptideSeq", "ProteinAccession")])

	#ver peptidos unicos en la muestra (peptidos con solo un hit en proteinAccession)
	tmp_uniques<-aggregate(ProteinAccession ~ PeptideSeq , data=tmp, FUN=paste, collapse=",")
	tmp_uniques$uniques<-str_count(paste(tmp_uniques$ProteinAccession),",")+1

	#merge con el objeto tmp para asignar el atributo "unico o no" a cada peptido (unico=1, no= >1)
	tmp_uniques$ProteinAccession<-NULL
	tmp<-merge(tmp,tmp_uniques, by="PeptideSeq")

	#genero columna count para hacer el spectral counting

	tmp$count <- paste(tmp$datfile, tmp$Query, sep = "_")
	tmp1 <- tmp[, c("count", "ProteinAccession")]
	tmp1<-unique(tmp1)
	#spectral counting

	results_sample_Counts <- summaryBy(ProteinAccession ~ count, data = tmp1, FUN = function(x) length(x))
	colnames(results_sample_Counts) <- c("count", "SpC")

	#merge con el objeto tmp para asignar columna SpC

	SPC_merged<-merge(tmp, results_sample_Counts, by="count")

	#me quedo solo con los peptidos unicos

	unicos<-SPC_merged[which(SPC_merged$uniques==1),]
	#cojo las longitudes de las proteinas
	results_sample_Counts<-merge(unicos,db_target_nAA, by.x="ProteinAccession", by.y="ProteinAccession")

	#elimino columnas no necesarias
	results_sample_Counts_end<-results_sample_Counts[,c("ProteinAccession","SpC","lengths")]

	results_sample_Counts$uNSAF <- (as.numeric(results_sample_Counts$SpC)/as.numeric(results_sample_Counts$lengths))/sum(as.numeric(results_sample_Counts$SpC)/as.numeric(results_sample_Counts$lengths))

	results_sample_Counts$sample <- sample
	if(log=="T"){
		results_sample_Counts$uNSAF<-log10(results_sample_Counts$uNSAF)
	}
	return(results_sample_Counts)

}

spectralCounting_Pepe_dNSAF<- function(results_sample, db_target_nAA,log){
	sample = paste(unique(results_sample$sample))
	cat("Calculating spectral counting ... ", sample," ... \n")
	tmp <- unique(results_sample[, c("datfile", "Query", "PeptideSeq","ProteinAccession")])
	##ver peptidos unicos en la muestra (peptidos con solo un hit en proteinAccession)

	tmp_uniques<-aggregate(ProteinAccession ~ PeptideSeq , data=tmp, FUN=paste, collapse=",")
	tmp_uniques$uniques<-str_count(paste(tmp_uniques$ProteinAccession),",")+1

	#merge con el objeto tmp para asignar el atributo "unico o compartido" a cada peptido (unico=1, compartido= >1)
	tmp_uniques$ProteinAccession<-NULL
	tmp<-merge(tmp,tmp_uniques, by="PeptideSeq")

	########## ver nº peptidos por proteina y asignar al objeto las longitudes de proteina
	tmp$count <- paste(tmp$datfile, tmp$Query, sep = "_")

	one_peptide_or_more<-tmp[,c("PeptideSeq", "ProteinAccession")]
	one_or_more_res<-summaryBy(PeptideSeq ~ ProteinAccession, data=one_peptide_or_more, FUN= function(x) length(x))
	names(one_or_more_res)<-c("ProteinAccession","n_peptides")
	tmp2<-tmp
	ALL_merged<-merge(tmp2,one_or_more_res, by="ProteinAccession", all.x=T)
	ALL_merged<-merge(ALL_merged,db_target_nAA, by="ProteinAccession")

	#genero columna count para hacer el spectral counting

	tmp1 <- tmp[, c("count", "ProteinAccession")]
	tmp1<-unique(tmp1)
	#spectral counting

	results_sample_Counts <- summaryBy(ProteinAccession ~ count, data = tmp1, FUN = function(x) length(x))
	colnames(results_sample_Counts) <- c("count", "SpC")

	#merge con el objeto ALL_merged para asignar columna SpC

	SPC_merged<-merge(ALL_merged, results_sample_Counts, by="count")

	# sum(uSpC totales) (constante para todo el dataset, para calcular factor d)
	ALL_merged_uniques<-SPC_merged[which(SPC_merged$uniques==1),]

	SpC_sum<-sum(ALL_merged_uniques$SpC)

 	SPC_merged$uniques[SPC_merged$uniques >1]<-0     #1 para los unicos, 0 para los shared
	Proteinlist<-unique(paste0(SPC_merged$ProteinAccession))
	subsets_output<-data.frame()
	for(i in 1:length(Proteinlist)){
		subset<-SPC_merged[which(SPC_merged$ProteinAccession==Proteinlist[i]),]
		if (nrow(subset)==1){
			if(subset$uniques==1){
				subset$dNSAF<-((as.numeric(subset$SpC)/as.numeric(SpC_sum))+ as.numeric(subset$SpC))/as.numeric(subset$lengths)
			}else{
				subset$dNSAF<- 0
			}
		}else if (nrow(subset)>1){
			if (sum(subset$uniques)==0){
				subset$dNSAF<- 0
			}else {
				For_SpC_shared<-subset[which(subset$uniques==0),]
				For_SpC_uniques<-subset[which(subset$uniques==1),]
				list_shared<-c()
				for (i in 1:nrow(For_SpC_shared)){
					if(nrow(For_SpC_shared) !=0){
					Spc<-as.numeric(paste(For_SpC_shared$SpC[i]))*(sum(For_SpC_uniques$SpC)/SpC_sum)+sum(For_SpC_uniques$SpC)
					list_shared<-c(list_shared,Spc)
					}else{
					Spc<-sum(For_SpC_uniques$SpC)
					list_shared<-c(list_shared, Spc)
					}
				}

				subset$dNSAF<-sum(list_shared)/as.numeric(paste(subset$lengths[1]))
			}
		}
		subsets_output<-rbind(subsets_output, subset)
	}
	subset_output_corrected<-subsets_output[,c("ProteinAccession", "dNSAF")]
	subset_output_corrected$sample<-sample
	if(log=="T"){
		subset_output_corrected$dNSAF<-log10(subset_output_corrected$dNSAF)
	}

	return(subset_output_corrected)
}



### create a results directory

resultDir <- paste0(experiment_Path, "/Results/")
y <- paste0("mkdir -p ", resultDir)
system(y)
y <- NULL

for(i in 1:length(file_path)){
	cat(file_path[i],"\n")
	results_sample<-get(load(file_path[i]))
	if(nrow(results_sample)==0){
		cat("empty file","\n")
		next
	}else{
		if(method=="NSAF"){
			NSAF<-spectralCounting_Pepe_NSAF(results_sample,db_target_nAA,log)
			namefile<-basename(file_path[i])
			new_name_F<-substr(paste(namefile),1,nchar(namefile)-4)
			results_file<-paste0(resultDir,new_name_F)
			save(NSAF, file=paste0(results_file,"NSAF_SC.rda"))   #en el analisis no tienen este nombre, se llaman igual que ei input
			NSAF$ID<-lapply(strsplit(paste(NSAF$ProteinAccession),"\\|"), "[", 2)
			write.table(paste(NSAF$ID), file=paste0(results_file,"_id.txt"), row.names=F, col.names=F, quote=F, sep="\t")
		}else if (method=="uNSAF"){
			uNSAF<-spectralCounting_Pepe_uNSAF(results_sample,db_target_nAA, log)
			namefile<-basename(file_path[i])
			new_name_F<-substr(paste(namefile),1,nchar(namefile)-4)
			results_file<-paste0(resultDir,new_name_F)
			#new_name<-paste0(new_name_F,"_spectral_counting_",method,".rda")
			save(uNSAF, file=paste0(results_file,"uNSAF_SC.rda"))

		}else if (method=="dNSAF"){
			dNSAF<-spectralCounting_Pepe_dNSAF(results_sample,db_target_nAA, log)
			namefile<-basename(file_path[i])
			new_name_F<-substr(paste(namefile),1,nchar(namefile)-4)
			#new_name<-paste0(new_name_F,"_spectral_counting_",method,".rda")
			results_file<-paste0(resultDir,new_name_F)
			save(dNSAF, file=paste0(results_file,"dNSAF_SC.rda"))

		}
	}
}

# selected_PROTEIN<-paste(uNSAF$ProteinAccession)
# NSAF_Selected<-NSAF[which(NSAF$ProteinAccession %in% selected_PROTEIN),]
# dNSAF_Selected<-dNSAF[which(dNSAF$ProteinAccession %in% selected_PROTEIN),]
# dNSAF_Selected$dNSAF[dNSAF_Selected$dNSAF ==-Inf]<-0
#
# dataframe_1<-merge(NSAF,dNSAF, by="ProteinAccession")
# dataframe_1<-merge(dataframe_1,uNSAF, by="ProteinAccession")
# final_df2<-dataframe_1[,c("ProteinAccession","NSAF","dNSAF")]
#
# pdf("/home/margaret/data/pepe/R_spectral_counting.pdf")
# upper.panel<-function(x, y){
#   points(x,y, pch=19, col=c("red", "green3", "blue")[final_df2$ProteinAccession])
#   r <- round(cor(x, y), digits=2)
#   txt <- paste0("R = ", r)
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   text(0.5, 0.9, txt)
# }
#
#
# pairs(~ NSAF + dNSAF, data=final_df2,
#       upper.panel = upper.panel)
# dev.off()
