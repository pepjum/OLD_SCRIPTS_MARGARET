##### Example   Rscript NEW_CLASSIFIER_SEP19_2.R /home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Dat_files/ PXD001383 mascot 5

args=(commandArgs(TRUE))

library(protr)
library(caret)
library(ggplot2)    # features dependientes de la secuencia del peptido
library(gridExtra)
library(stringr)
library(plyr)
library(dplyr)

source("/home/margaret/data/pepe/scripts/functions_Classifier.R")

path<-args[1] #"/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Dat_files/"
dataset<-args[2]   #"PXD001383"
SE<-args[3]  #"mascot"
iteration_n<-args[4]
#decoy_id<-args[5]    #"DECOY"
#minn_AA <- args[6]   # 9
#maxn_AA <- args[7]  # 30

psmFDRvalue = 0.01
protFDRvalue = 0.01
pepScore="score"

prot_id = "database"
pep_id="Pep_scan"
pep_col="PeptideSeq"

if(SE=="mascot"){

    files_target_n<-list.files(paste0(path,dataset), pattern = ".dat")
    files_decoy_n<-list.files(paste0(path,paste0(dataset,"-D/")),pattern=".dat")

}else if(SE=="xtandem"){

    files_target_n<-list.files(paste0(path,dataset), pattern = "_corrected.tsv")
    files_decoy_n<-list.files(paste0(path,paste0(dataset,"-D/")),pattern="_corrected.tsv")

}else if(SE=="comet"){

    files_target_n<-list.files(paste0(path,dataset), pattern = "_corrected.txt")
    files_decoy_n<-list.files(paste0(path,paste0(dataset,"-D/")),pattern="_corrected.txt")

}

files_target<-paste(path,dataset,files_target_n,sep="/")
experiment_target<-data.frame("V1"=paste(files_target),"V2"=paste("TARGET"),"V3"=paste(dataset),"V4"=seq(1:length(files_target)),"V5"=paste(files_target_n))
files_decoy<-paste(path,paste0(dataset,"-D"),files_decoy_n,sep="/")
experiment_decoy<-data.frame("V1"=paste(files_decoy),"V2"=paste("DECOY"),"V3"=paste0(dataset,"-D"),"V4"=seq(1:length(files_decoy)),"V5"=paste(files_decoy_n))
if(nrow(experiment_decoy)==nrow(experiment_target)){

    experimentData<-rbind(experiment_target,experiment_decoy)

}else{

    cat("WARNING: Number of fractions in Target folder is not equal to Decoy Folder. You may be interested in reviewing searches","\n")
    cat("Continuing the analysis","\n")
    experimentData<-rbind(experiment_target,experiment_decoy)
}



cat("Loading files \n")

dataPSMMat <- data.frame()
j=1

for (j in 1:nrow(experimentData))

{
	if (j == 1)
	{
		if (SE == "mascot")
		{
			dataPSMMat <-  data.frame("file" = paste(experimentData[1,1]), parseMascotStandardOutput(paste(experimentData[1,1])), "database" = paste(experimentData[1,2]), "sample" = paste(experimentData[1,3]), "fraction" = paste(experimentData[1,4]), "datfile" = paste(experimentData[1,5]))
		}
		else if (SE == "comet")
		{
			dataPSMMat <-  data.frame("file" = paste(experimentData[1,1]), parseCometStandardOutput(paste(experimentData[1,1])), "database" = paste(experimentData[1,2]), "sample" = paste(experimentData[1,3]), "fraction" = paste(experimentData[1,4]), "datfile" = paste(experimentData[1,5]))
		}
		else if (SE == "xtandem")
		{
			dataPSMMat <-  data.frame("file" = paste(experimentData[1,1]), parseXTandemStandardOutput(paste(experimentData[1,1])), "database" = paste(experimentData[1,2]), "sample" = paste(experimentData[1,3]), "fraction" = paste(experimentData[1,4]), "datfile" = paste(experimentData[1,5]))
		}

		else
		{
			cat("Error. Specify search engine parameter correctly. mascot for Mascot. comet for Comet. tandem for X!tandem")
			break
		}
	}else{
		if (SE == "mascot")
		{
			dataPSMMat <- rbind(dataPSMMat, data.frame("file" = paste(experimentData[j,1]), parseMascotStandardOutput(paste(experimentData[j,1])), "database" = paste(experimentData[j,2]), "sample" = paste(experimentData[j,3]), "fraction" = paste(experimentData[j,4]), "datfile" = paste(experimentData[j,5])))
		}
		else if (SE == "comet")
		{
			dataPSMMat <- rbind(dataPSMMat, data.frame("file" = paste(experimentData[j,1]), parseCometStandardOutput(paste(experimentData[j,1])), "database" = paste(experimentData[j,2]), "sample" = paste(experimentData[j,3]), "fraction" = paste(experimentData[j,4]), "datfile" = paste(experimentData[j,5])))
		}
		else if (SE == "xtandem")
		{
			dataPSMMat <- rbind(dataPSMMat, data.frame("file" = paste(experimentData[j,1]), parseXTandemStandardOutput(paste(experimentData[j,1])), "database" = paste(experimentData[j,2]), "sample" = paste(experimentData[j,3]), "fraction" = paste(experimentData[j,4]), "datfile" = paste(experimentData[j,5])))
		}
	}
}

cat("computing some features...\n")

dataPSMMat$PepLen<-nchar(as.character(dataPSMMat$PeptideSeq))

dataPSMMat$database<-as.factor(dataPSMMat$database)

for(number in c(4:6,9,19)){
    dataPSMMat[,number]<-as.numeric(dataPSMMat[,number])
}

cat("preparing training set....\n")

ALL_DECOY<-ALL_DECOY<-dataPSMMat[which(dataPSMMat$database=="DECOY"),]

DECOYS_for_TRAINING<-sample_n(ALL_DECOY, nrow(ALL_DECOY)/2)

ALL_TARGET<-dataPSMMat[which(dataPSMMat$database=="TARGET"),]
TMP_OBJECT<-rbind(ALL_TARGET, DECOYS_for_TRAINING)    #ESTE objecto es sobre el que hay que hacer las iteraciones

DF_with_FDR<-psmFDR(TMP_OBJECT, "score", "DECOY",concat_decoy=0, prot_id="database")

DF_FDR_PASSED<-DF_with_FDR[which(DF_with_FDR$psmFDR<=0.01),]


if(nrow(DF_FDR_PASSED==0)){

    cat("There are some DECOY PSMs hits with scores higher than the best TARGET PSM. Could not create training dataset. Please remove these hits and try again...\n")
    break
}

DF_FDR_PASSED$psmFDR<-NULL

df_selected<-rbind(DF_FDR_PASSED,DECOYS_for_TRAINING)   #con este entrenamos la primera glm
trainDescr<-df_selected

cat("plotting some characteristics before training classificator...\n")

list_plots<-list()
counts<-1
TARGET<-trainDescr[which(trainDescr$database=="TARGET"),]
DECOY<-trainDescr[which(trainDescr$database=="DECOY"),]

for (i in c(5,6,9,19)){
    plotter<-data.frame()
    plotter<-rbind(plotter, TARGET[,c(15,i)])
    plotter<-rbind(plotter, DECOY[,c(15,i)])
    plotter[,2]<-as.numeric(plotter[,2])
    plotter_aes<-names(plotter[2])

    plot<-ggplot(plotter, aes_string(x="database", y=plotter_aes, color="database"))+geom_violin()+ggtitle(plotter_aes)

    list_plots[[counts]]<-plot
    counts<-counts+1
}

ggsave(paste0(arg1,"BEFORE_CLASSIFICATIONS_DIFFERENCIES.pdf"), arrangeGrob(grobs = list_plots))

set.seed(1)
ctrl <- trainControl(method="repeatedcv", repeats=3, classProbs=TRUE, summaryFunction= twoClassSummary)

psms<-c()
peptides<-c()
proteins<-c()
dataframes<-list()
glms<-list()

psms_all<-c()
peptides_all<-c()
proteins_all<-c()
dataframes_all<-list()
glms_all<-list()

psms_all_clean<-c()
peptides_all_clean<-c()
proteins_all_clean<-c()


counter<-1

#numero de iteraciones
trainDescr<-trainDescr[,c(1,2,3,4,15,5,6,7,8,9,10,11,12,13,14,16,17,18,19)]
TMP_OBJECT<-TMP_OBJECT[,c(1,2,3,4,15,5,6,7,8,9,10,11,12,13,14,16,17,18,19)]

for (i in seq(1:iteration_n)){

    if (i==1){
        cat("training method iteration ", i,"\n")
        glmFit<-train(database~., data=trainDescr[,c(5,6,7,10,19)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

        cat("reranking test data", i,"\n")

        iteration<-predict(glmFit, TMP_OBJECT[,c(5,6,7,10,19)], type = 'prob')  #primera
        TMP_OBJECT[,20]<-as.numeric(iteration[,2])
        names(TMP_OBJECT)[20]<-paste0("p_TARGET_",i)
        a<-names(TMP_OBJECT[20])
        FIRST_ITERATION_FDR<-psmFDR(TMP_OBJECT, a, "DECOY",concat_decoy=0, prot_id="Proteins")
        df_temp<-FIRST_ITERATION_FDR_passed[,c(1:19)]

        FIRST_ITERATION_FDR_passed<-FIRST_ITERATION_FDR[which(FIRST_ITERATION_FDR$psmFDR <=0.01),]
        dataframes[[counter]]<-df_temp
        psms[counter]<-nrow(FIRST_ITERATION_FDR_passed)
        peptides[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed$Peptide)))
        proteins[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed$Proteins)))
        glms[[counter]]<-glmFit

        cat("applying to full data ",i,"\n")

        iteration_all<-predict(glmFit, dataPSMMat[,c(5,6,7,10,19)], type = 'prob')  #primera
        df[,20]<-as.numeric(iteration_all[,2])
        names(df)[20]<-paste0("p_TARGET_",i)
        a<-names(df[20])
        FIRST_ITERATION_FDR_all<-psmFDR(df, a, "DECOY",concat_decoy=0, prot_id="database")
        FIRST_ITERATION_FDR_passed_all<-FIRST_ITERATION_FDR_all[which(FIRST_ITERATION_FDR_all$psmFDR <=0.01),]
        df_temp_all<-FIRST_ITERATION_FDR_passed_all[,c(1:19)]

        dataframes_all[[counter]]<-df_temp_all
        psms_all[counter]<-nrow(FIRST_ITERATION_FDR_passed_all)
        peptides_all[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed_all$Peptide)))
        proteins_all[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed_all$Proteins)))

        only_targets<-FIRST_ITERATION_FDR_passed_all[which(FIRST_ITERATION_FDR_passed_all$Label=="TARGET"),]
        psms_all_clean[counter]<-nrow(only_targets)
        peptides_all_clean[counter]<-length(unique(paste(only_targets$Peptide)))
        proteins_all_clean[counter]<-length(unique(paste(only_targets$Proteins)))

        glms_all[[counter]]<-glmFit

        counter<-counter+1

    }else{

        trainDescr_tmp<-rbind(dataframes[[counter-1]], DECOYS_for_TRAINING)

        tmp_train_iter_FDR<-psmFDR(trainDescr_tmp, "score", "DECOY",concat_decoy=0, prot_id="database")   #creo que no es necesario

        tmp_train_iter_FDR_passed<-tmp_train_iter_FDR[which(tmp_train_iter_FDR$psmFDR<=0.01),]   # creo q  no es necesario
        tmp_train_iter_FDR_passed$psmFDR<-NULL

        train_DESCR_2<-rbind(tmp_train_iter_FDR_passed,DECOYS_for_TRAINING)

        cat("training method iteration ", i,"\n")

        glmFit_2<-train(database~., data=train_DESCR_2[,c(5,6,7,10,19)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

        glms[[counter]]<-glmFit_2
        glms_all[[counter]]<-glmFit_2

        cat("reranking test data", i,"\n")

        iteration<-predict(glmFit_2, TMP_OBJECT[,c(5,6,7,10,19)], type = 'prob')
        TMP_OBJECT[,19+i]<-as.numeric(iteration[,2])
        names(TMP_OBJECT)[19+i]<-paste0("p_TARGET_",i)
        a<-names(TMP_OBJECT)[19+i]
        ITERATION_FDR<-psmFDR(TMP_OBJECT, a, "DECOY",concat_decoy=0, prot_id="database")
        ITERATION_FDR_passed<-ITERATION_FDR[which(FIRST_ITERATION_FDR$psmFDR <=0.01),]
        df_temp<-ITERATION_FDR_passed[,c(1:19)]

        dataframes[[counter]]<-df_temp
        psms[counter]<-nrow(df_temp)
        peptides[counter]<-length(unique(paste(ITERATION_FDR_passed$Peptide)))
        proteins[counter]<-length(unique(paste(ITERATION_FDR_passed$Proteins)))

        cat("applying to full data", i,"\n")

        iteration_all<-predict(glmFit_2, dataPSMMat[,c(5,6,7,10,19)], type = 'prob')
        df[,19+i]<-as.numeric(iteration_all[,2])
        names(df)[19+i]<-paste0("p_TARGET_",i)
        a<-names(df)[19+i]
        ITERATION_FDR_all<-psmFDR(df, a, "DECOY",concat_decoy=0, prot_id="database")
        ITERATION_FDR_passed_all<-ITERATION_FDR_all[which(FIRST_ITERATION_FDR_all$psmFDR <=0.01),]
        df_temp_all<-ITERATION_FDR_passed_all[,c(1:19)]

        dataframes_all[[counter]]<-df_temp_all
        psms_all[counter]<-nrow(df_temp_all)
        peptides_all[counter]<-length(unique(paste(ITERATION_FDR_passed_all$Peptide)))
        proteins_all[counter]<-length(unique(paste(ITERATION_FDR_passed_all$Proteins)))

        only_targets<-df_temp_all[which(df_temp_all$Label=="TARGET"),]
        psms_all_clean[counter]<-nrow(only_targets)
        peptides_all_clean[counter]<-length(unique(paste(only_targets$Peptide)))
        proteins_all_clean[counter]<-length(unique(paste(only_targets$Proteins)))


        counter<-counter+1
    }
}

cat("computing psmFDR at interation 0..\n")   #en el FULL OBJECT

df_FDR_iteration_0<-psmFDR(df, pepScore="score", decoy_id="DECOY",concat_decoy=0, prot_id="Label")

df_FDR_iteration_0_passed<-df_FDR_iteration_0[which(df_FDR_iteration_0$psmFDR <=0.01),]

tmp_vect_psm_0<-nrow(df_FDR_iteration_0_passed)
tmp_peptides_0<-length(unique(paste(df_FDR_iteration_0_passed$Peptide)))
tmp_proteins_0<-length(unique(paste(df_FDR_iteration_0_passed$Proteins)))

df_FDR_iteration_0_passed_clean<-df_FDR_iteration_0_passed[which(df_FDR_iteration_0_passed$Label=="TARGET"),]

tmp_vect_psm_0_clean<-nrow(df_FDR_iteration_0_passed_clean)
tmp_peptides_0_clean<-length(unique(paste(df_FDR_iteration_0_passed_clean$Peptide)))
tmp_proteins_0_clean<-length(unique(paste(df_FDR_iteration_0_passed_clean$Proteins)))

#totales con decoys que pasan fdr
total_psms<-append(tmp_vect_psm_0,psms_all)
total_peptides<-append(tmp_peptides_0,peptides_all)
total_proteins<-append(tmp_proteins_0,proteins_all)

### totales solo targets

total_psms_clean<-append(tmp_vect_psm_0_clean,psms_all_clean)
total_peptides_clean<-append(tmp_peptides_0_clean,peptides_all_clean)
total_proteins_clean<-append(tmp_proteins_0_clean,proteins_all_clean)

cat("plotting results of the improvements in the FDR calculation.....\n")

plotter<-data.frame("iteration"=seq(0:iteration_n)-1,"psms"=paste(total_psms_clean), "peptides"=paste(total_peptides_clean),"proteins"=paste(total_proteins_clean))

list_plots<-list()
counts<-1
for (i in (2:4)){
    plotter_tmp<-data.frame()
    plotter_tmp<-rbind(plotter_tmp, plotter[,c(1,i)])
    plotter_aes<-names(plotter_tmp[2])

    plot<-ggplot(plotter_tmp, aes_string(x="iteration", y=plotter_aes, group=1))+geom_point() +geom_line() +ggtitle(plotter_aes) + coord_cartesian(xlim = c(0, 5))

    list_plots[[counts]]<-plot
    counts<-counts+1
}


cat("end of the script...\n")
