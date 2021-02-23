library(protr)
library(caret)
library(ggplot2)    # features dependientes de la secuencia del peptido
library(gridExtra)
library(stringr)
library(plyr)
library(dplyr)

psmFDR <- function(protPep, pepScore="hyperscore", decoy_id="DECOY", concat_decoy=0, prot_id="Label")
{
	protPep_decoy <- protPep[order(as.numeric(paste(protPep[,pepScore])), decreasing=TRUE),]
	decoyVec <- vector(mode = "numeric", nrow(protPep_decoy))
	decoyVec[grep(decoy_id, protPep_decoy[,prot_id])] <- 1
	if(concat_decoy==1)
	{
		protPep_decoy$psmFDR <- (2*cumsum(decoyVec))/(1:nrow(protPep_decoy))
	}else
	{
		protPep_decoy$psmFDR <- cumsum(decoyVec)/((1:nrow(protPep_decoy))-(cumsum(decoyVec)))
	}
	# todo lo mayor a 1, lo ponemos a 1, ya que no debería poder darse
	protPep_decoy[protPep_decoy$psmFDR > 1, "psmFDR"] <- 1
	return(protPep_decoy)

}

if (SE=="tandem"){
    df<-read.table("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Tandem_Files/PXD001383/percolator/ALL_XTANDEM.pin", header=T, fill=NA,stringsAsFactors=F)
    df<-read.table("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Tandem_Files/PXD001468/percolator/ALL_XTANDEM.pin", header=T, fill=NA,stringsAsFactors=F)
    
    df[which(df$Proteins==""),]<-NA
    df$Label<-as.numeric(df$Label)
    df<-df[complete.cases(df),]
    df$Label[which(df$Label==-1)]<-"DECOY"
    df$Label[which(df$Label==1)]<-"TARGET"
}else if (SE=="mascot"){
    df<-get(load("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Dat_files/results_Peptides_mascot_PXD001383_dataPSMMat.rda"))
    df$Label<-NA
    df$Label[which(df$database=="T")]<-"TARGET"
    df$Label[which(df$database=="D")]<-"DECOY"
}
#preparar los datos


if(SE=="tandem"){

    cat("tandem","\n")

    for (i in (4:17)){

        df[,i]<-as.numeric(df[,i])

    }

    df[,2]<-as.factor(df[,2])

}else if(SE=="mascot"){

    cat("mascot","\n")

    for (number in c(4,5,6,9)){

        df[,number]<-as.numeric(df[,number])

    }

    df$Label<-as.factor(df$Label)

}



#### DIVISION DE LAS DECOYS EN DOS SUBSETS (con una seleccion aleatoria)

ALL_DECOY<-ALL_DECOY<-df[which(df$Label=="DECOY"),]
DECOYS_for_TRAINING<-sample_n(ALL_DECOY, nrow(ALL_DECOY)/2)
ALL_TARGET<-df[which(df$Label=="TARGET"),]

#crear objeto con targets y la mitad de DECOYS

TMP_OBJECT<-rbind(ALL_TARGET, DECOYS_for_TRAINING)    #ESTE objecto es sobre el que hay que hacer las iteraciones

if(SE=="tandem"){

    DF_with_FDR<-psmFDR(TMP_OBJECT, pepScore="hyperscore", decoy_id="DECOY",concat_decoy=0, prot_id="Label")

}else if(SE=="mascot"){

    DF_with_FDR<-psmFDR(TMP_OBJECT, pepScore="score", decoy_id="DECOY",concat_decoy=0, prot_id="Label")

}
#selecciono solo los que pasan el FDR
DF_FDR_PASSED<-DF_with_FDR[which(DF_with_FDR$psmFDR<=0.01),]
DF_FDR_PASSED$psmFDR<-NULL

#añado de nuevo el mismo set de decoys a lo que pasa el FDR para tener BEST_HITS
df_selected<-rbind(DF_FDR_PASSED,DECOYS_for_TRAINING)   #con este entrenamos la primera glm
trainDescr<-df_selected

#control para el clasificador y definicion de variables de salida
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
iteration_n<-5

for (i in seq(1:iteration_n)){

    if (i==1){

        cat("training method iteration ", i,"\n")

        if(SE=="tandem"){

            glmFit<-train(Label~., data=trainDescr[,c(2,4,5,7)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

            cat("applying to test data", i,"\n")

            iteration<-predict(glmFit, TMP_OBJECT[,c(2,4,5,7)], type = 'prob')  #primera
            TMP_OBJECT[,20]<-as.numeric(iteration[,2])
            names(TMP_OBJECT)[20]<-paste0("p_TARGET_",i)
            a<-names(TMP_OBJECT[20])
            FIRST_ITERATION_FDR<-psmFDR(TMP_OBJECT, a, "DECOY",concat_decoy=0, prot_id="Proteins")
            df_temp<-FIRST_ITERATION_FDR_passed[,c(1:19)]

        }else if(SE=="mascot"){

            glmFit<-train(Label~., data=trainDescr[,c(4,5,6,9,20)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

            cat("applying to test data", i,"\n")

            iteration<-predict(glmFit, TMP_OBJECT[,c(4,5,6,9,20)], type = 'prob')  #primera
            TMP_OBJECT[,21]<-as.numeric(iteration[,2])
            a<-names(TMP_OBJECT)[21]
            FIRST_ITERATION_FDR<-psmFDR(TMP_OBJECT, a, "DECOY",concat_decoy=0, prot_id="Proteins")
            df_temp<-FIRST_ITERATION_FDR_passed[,c(1:20)]

        }

        FIRST_ITERATION_FDR_passed<-FIRST_ITERATION_FDR[which(FIRST_ITERATION_FDR$psmFDR <=0.01),]
        dataframes[[counter]]<-df_temp
        psms[counter]<-nrow(FIRST_ITERATION_FDR_passed)
        peptides[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed$Peptide)))
        proteins[counter]<-length(unique(paste(FIRST_ITERATION_FDR_passed$Proteins)))
        glms[[counter]]<-glmFit

        if(SE=="tandem"){

            cat("applying to full data ",i,"\n")

            iteration_all<-predict(glmFit, df[,c(2,4,5,7)], type = 'prob')  #primera
            df[,20]<-as.numeric(iteration_all[,2])
            names(df)[20]<-paste0("p_TARGET_",i)
            a<-names(df[20])
            FIRST_ITERATION_FDR_all<-psmFDR(df, a, "DECOY",concat_decoy=0, prot_id="Proteins")
            FIRST_ITERATION_FDR_passed_all<-FIRST_ITERATION_FDR_all[which(FIRST_ITERATION_FDR_all$psmFDR <=0.01),]
            df_temp_all<-FIRST_ITERATION_FDR_passed_all[,c(1:19)]

        }else if(SE=="mascot"){

            cat("applying to full data ",i,"\n")

            iteration_all<-predict(glmFit, df[,c(4,5,6,9,20)], type = 'prob')  #primera
            df[,21]<-as.numeric(iteration_all[,2])
            names(df)[21]<-paste0("p_TARGET_",i)
            a<-names(df[21])
            FIRST_ITERATION_FDR_all<-psmFDR(df, a, "DECOY",concat_decoy=0, prot_id="Proteins")
            FIRST_ITERATION_FDR_passed_all<-FIRST_ITERATION_FDR_all[which(FIRST_ITERATION_FDR_all$psmFDR <=0.01),]
            df_temp_all<-FIRST_ITERATION_FDR_passed_all[,c(1:20)]

        }

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

        if(SE=="tandem"){

            tmp_train_iter_FDR<-psmFDR(trainDescr_tmp, "hyperscore", "DECOY",concat_decoy=0, prot_id="Label")

        }else if(SE=="mascot"){

            tmp_train_iter_FDR<-psmFDR(trainDescr_tmp, "score", "DECOY",concat_decoy=0, prot_id="Label")

        }

        tmp_train_iter_FDR_passed<-tmp_train_iter_FDR[which(tmp_train_iter_FDR$psmFDR<=0.01),]
        tmp_train_iter_FDR_passed$psmFDR<-NULL

        train_DESCR_2<-rbind(tmp_train_iter_FDR_passed,DECOYS_for_TRAINING)

        cat("training method iteration ", i,"\n")

        if(SE=="tandem"){

            glmFit_2<-train(Label~., data=train_DESCR_2[,c(2,4,5,7)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

        }else if(SE=="mascot"){

            glmFit_2<-train(Label~., data=train_DESCR_2[,c(4,5,6,9,20)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

        }

        glms[[counter]]<-glmFit_2
        glms_all[[counter]]<-glmFit_2

        cat("applying to test data", i,"\n")

        if(SE=="tandem"){

            iteration<-predict(glmFit_2, TMP_OBJECT[,c(2,4,5,7)], type = 'prob')
            TMP_OBJECT[,19+i]<-as.numeric(iteration[,2])
            names(TMP_OBJECT)[19+i]<-paste0("p_TARGET_",i)
            a<-names(TMP_OBJECT)[19+i]
            ITERATION_FDR<-psmFDR(TMP_OBJECT, a, "DECOY",concat_decoy=0, prot_id="Proteins")
            ITERATION_FDR_passed<-ITERATION_FDR[which(FIRST_ITERATION_FDR$psmFDR <=0.01),]
            df_temp<-ITERATION_FDR_passed[,c(1:19)]

        }else if(SE=="mascot"){

            iteration<-predict(glmFit_2, TMP_OBJECT[,c(4,5,6,9,20)], type = 'prob')
            TMP_OBJECT[,20+i]<-as.numeric(iteration[,2])
            names(TMP_OBJECT)[20+i]<-paste0("p_TARGET_",i)
            a<-names(TMP_OBJECT)[20+i]
            ITERATION_FDR<-psmFDR(TMP_OBJECT, a, "DECOY",concat_decoy=0, prot_id="Proteins")
            ITERATION_FDR_passed<-ITERATION_FDR[which(FIRST_ITERATION_FDR$psmFDR <=0.01),]
            df_temp<-ITERATION_FDR_passed[,c(1:20)]

        }

        dataframes[[counter]]<-df_temp
        psms[counter]<-nrow(df_temp)
        peptides[counter]<-length(unique(paste(ITERATION_FDR_passed$Peptide)))
        proteins[counter]<-length(unique(paste(ITERATION_FDR_passed$Proteins)))

        cat("applying to full data", i,"\n")

        if(SE=="tandem"){

            iteration_all<-predict(glmFit_2, df[,c(2,4,5,7)], type = 'prob')
            df[,19+i]<-as.numeric(iteration_all[,2])
            names(df)[19+i]<-paste0("p_TARGET_",i)
            a<-names(df)[19+i]
            ITERATION_FDR_all<-psmFDR(df, a, "DECOY",concat_decoy=0, prot_id="Proteins")
            ITERATION_FDR_passed_all<-ITERATION_FDR_all[which(FIRST_ITERATION_FDR_all$psmFDR <=0.01),]
            df_temp_all<-ITERATION_FDR_passed_all[,c(1:19)]

        }else if(SE=="mascot"){

            iteration_all<-predict(glmFit_2, df[,c(4,5,6,9,20)], type = 'prob')
            df[,20+i]<-as.numeric(iteration_all[,2])
            names(df)[20+i]<-paste0("p_TARGET_",i)
            a<-names(df)[20+i]
            ITERATION_FDR_all<-psmFDR(df, a, "DECOY",concat_decoy=0, prot_id="Proteins")
            ITERATION_FDR_passed_all<-ITERATION_FDR_all[which(FIRST_ITERATION_FDR_all$psmFDR <=0.01),]
            df_temp_all<-ITERATION_FDR_passed_all[,c(1:20)]

        }

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

if(SE=="tandem"){

    df_FDR_iteration_0<-psmFDR(df, pepScore="hyperscore", decoy_id="DECOY",concat_decoy=0, prot_id="Label")

}else if(SE=="mascot"){

    df_FDR_iteration_0<-psmFDR(df, pepScore="score", decoy_id="DECOY",concat_decoy=0, prot_id="Label")

}

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

### GUARDAR TODAS LOS PLOTS
ggsave("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Results_classificator_psm_pep_prot_only_TARGETS.pdf", arrangeGrob(grobs = list_plots))

plotter<-data.frame("iteration"=seq(0:iteration_n)-1,"psms"=paste(total_psms), "peptides"=paste(total_peptides),"proteins"=paste(total_proteins))

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

### GUARDAR TODAS LOS PLOTS
ggsave("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Results_classificator_psm_pep_prot_all_FDR_passed.pdf", arrangeGrob(grobs = list_plots))
