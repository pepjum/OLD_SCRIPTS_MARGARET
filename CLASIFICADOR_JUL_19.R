###CLASIFICADORES
library(protr)
library(caret)
library(ggplot2)    # features dependientes de la secuencia del peptido
library(gridExtra)
library(stringr)
library(plyr)
library(dplyr)
# preparar datos de XTANDEM
df<-read.table("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/Tandem_Files/PXD001383/percolator/ALL_XTANDEM.pin", header=T, fill=NA,stringsAsFactors=F)
df[which(df$Proteins==""),]<-NA
df$Label<-as.numeric(df$Label)
df<-df[complete.cases(df),]

#preparar los datos

df$Label[which(df$Label==-1)]<-"DECOY"
df$Label[which(df$Label==1)]<-"TARGET"

for (i in (4:17)){
    df[,i]<-as.numeric(df[,i])
}
df[,2]<-as.factor(df[,2])


#### DIVISION DE LAS DECOYS EN DOS SUBSETS (con una seleccion aleatoria)

ALL_DECOY<-ALL_DECOY<-df[which(df$Label=="DECOY"),]

DECOYS_for_TRAINING<-sample_n(ALL_DECOY, nrow(ALL_DECOY)/2)

ALL_TARGET<-df[which(df$Label=="TARGET"),]
TMP_OBJECT<-rbind(ALL_TARGET, DECOYS_for_TRAINING)    #ESTE objecto es sobre el que hay que hacer las iteraciones

#TMP_OBJECT es el objeto con todo target y la mitad de decoys, a esto le calculo FDR, me quedo con lo menor al 1% y le vuelvo a añadir la mitad de las decoys

#CALCULAR EL FDR CON LA FUNCION DE ALBA

psmFDR <- function(protPep, pepScore="hyperscore", decoy_id="DECOY", concat_decoy=0, prot_id="prot_acc")
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


DF_with_FDR<-psmFDR(TMP_OBJECT, "score", "DECOY",concat_decoy=0, prot_id="database")

#selecciono solo los que pasan el FDR
DF_FDR_PASSED<-DF_with_FDR[which(DF_with_FDR$psmFDR<=0.01),]
#1918

DF_FDR_PASSED$psmFDR<-NULL

df_selected<-rbind(DF_FDR_PASSED,DECOYS_for_TRAINING)   #con este entrenamos la primera glm

################################## variables de PERCOLATOR INTERNAS

# df_selected$UniqID<-rownames(df_selected)
# # la diferencia entre un peptido target y uno decoy está en la secuencia. Tengo que transfomar la secuencia string en una lista numérica, por lo que tengo que crear un diccionario de valores para transformarlo:
# findProteinMatches<-function(peptide,df){
#         tmp<-df[which(df$Peptide==peptide),]
#         list<-unlist(tmp$Proteins)
#         list<-unique(list)
#     return(list)
# }
#
# intrasetFeatures<-function(df){
#     numberPeptides<-list()
#     uniqPep<-list()
#     numberProteins<-list()
#     for (row in 1:nrow(df)){
#         peptide<-paste(df[row, 'Peptide'])
#         numberPeptides[peptide]<-0
#
#         for(Protein in findProteinMatches(peptide, df)){
#             numberProteins[Protein]<-0
#             uniqPep[Protein]<-NULL
#         }
#     }
#     for (row in 1:nrow(df)){
#         peptide<-paste(df[row,"Peptide"])
#         numberPeptides[peptide]<-as.numeric(numberPeptides[peptide])+1
#         for(Protein in findProteinMatches(peptide,df)){
#             numberProteins[Protein]<-as.numeric(numberProteins[Protein])+1
#             if (!( peptide %in% uniqPep[Protein])){
#                 uniqPep[[Protein]]<-c(uniqPep[[Protein]],peptide)
#             }
#         }
#     }
#     exit<-data.frame(numPep=NULL, numProt=NULL, pepSite=NULL)
#     for (row in 1:nrow(df)){
#          peptide<-paste(df[row,"Peptide"])
#          numPep<-numberPeptides[peptide]
#          numProt<-0
#          pepSite<-0
#          for(Protein in findProteinMatches(peptide,df)){
#              numProt<-max(numProt, numberProteins[[Protein]])
#              pepSite<-max(pepSite,length(unique(uniqPep[[Protein]])))
#          }
#     exit<-rbind( exit, data.frame(numPep=log(numberPeptides[[peptide]]),numProt=log(numProt),pepSite=log(pepSite)))
#     }
#     exit$UniqID <- df$UniqID
#     return(exit)
# }
# output<-intrasetFeatures(df_selected)
#
############################################### LO ANTERIOR NO LO HAGO EN ESTA PRUEBA


# df_selected<-df_selected[,c(1:2,4,5,7,18)]

#object_to_train<-merge(df_selected, output, by="UniqID")
object_to_train<-df_selected


### Voy a hacer boxplots de las caracteristicas que he seleccionado para ver si hay diferencias significativas entre T y D

TARGET<-object_to_train[which(object_to_train$Label=="TARGET"),]
DECOY<-object_to_train[which(object_to_train$Label=="DECOY"),]
TARGET$Peptide<-NULL
DECOY$Peptide<-NULL
#
#
# #ejemplos de t.test
# t.test(as.numeric(TARGET$PepLen), as.numeric(DECOY$PepLen))$p.value
# t.test(as.numeric(TARGET$numPep), as.numeric(DECOY$numPep))$p.value
#
#
##########
#PLOTS DE LAS CARACTERISTICAS
###################

list_plots<-list()
counts<-1
for (i in (4:7)){
    plotter<-data.frame()
    plotter<-rbind(plotter, TARGET[,c(2,i)])
    #plotter$Label<-"TARGET"
    plotter<-rbind(plotter, DECOY[,c(2,i)])
    #plotter$Label[which(plotter$Label==-1)]<-"DECOY"
    plotter[,2]<-as.numeric(plotter[,2])
    plotter_aes<-names(plotter[2])

    plot<-ggplot(plotter, aes_string(x="Label", y=plotter_aes, color="Label"))+geom_violin()+ggtitle(plotter_aes)

    list_plots[[counts]]<-plot
    counts<-counts+1
}


### GUARDAR TODAS LOS PLOTS
ggsave("/home/margaret/data/pepe/16_CLASIFICADORES_JUL19/FEATURES_SELECTED_TARGET_VS_DECOY_PERCOLATOR_FDR.pdf", arrangeGrob(grobs = list_plots))

set.seed(1)

# object_to_train$Label[which(object_to_train$Label==-1)]<-"DECOY"
# object_to_train$Label[which(object_to_train$Label==1)]<-"TARGET"

# object_to_train$Label<-as.factor(object_to_train$Label)
#
# for (i in 3:5){
#     object_to_train[,i]<-as.numeric(object_to_train[,i])
# }
#

#preparar columnas de trainDescr

trainDescr<-object_to_train
# for (i in 3:17){
#     trainDescr[,i]<-as.numeric(trainDescr[,i])
# }


##### BUILDING MODEL 1ª iteracion

ctrl <- trainControl(method="repeatedcv", repeats=3, classProbs=TRUE, summaryFunction= twoClassSummary)
glmFit<-train(Label~., data=trainDescr[,c(2,4,5,7)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

#Clasificadores de ELI son: svmRadial, PLS, RPART, CS, GLM, JRIP, RF, NB, NNET, PART

############################################ 1º iteracion: clasificar el dataset original ################################################

first_iteration<-predict(glmFit, TMP_OBJECT[,c(2,4,5,7)], type = 'prob')  #primera

TMP_OBJECT$p_Target_first<-as.numeric(first_iteration[,2])

##### calcular el PSMFDR, ordenando en este caso por la probabilidad que nos ha dado el clasificador

FIRST_ITERATION_FDR<-psmFDR(TMP_OBJECT, "p_Target_first", "DECOY",concat_decoy=0, prot_id="Proteins")

FIRST_ITERATION_FDR_passed<-FIRST_ITERATION_FDR[which(FIRST_ITERATION_FDR$psmFDR <=0.01),]

FIRST_ITERATION_FDR_passed$psmFDR<-NULL

######################################## 2º iteracion: CREAR nueva glm. El objeto de partida será el de la primera iteracion que ha pasado el fdr + DECOYS_fOR_TRAINING, rerankeado por FDR  #########################################################################
FIRST_ITERATION_FDR_passed$p_Target_first<-NULL
trainDESCR_2_tmp<-rbind(FIRST_ITERATION_FDR_passed, DECOYS_for_TRAINING)

tmp_train_second_iter_FDR<-psmFDR(trainDESCR_2_tmp, "hyperscore", "DECOY",concat_decoy=0, prot_id="Proteins")
tmp_train_second_iter_FDR_passed<-tmp_train_second_iter_FDR[which(tmp_train_second_iter_FDR$psmFDR<=0.01),]
tmp_train_second_iter_FDR_passed$psmFDR<-NULL

train_DESCR_2<-rbind(tmp_train_second_iter_FDR_passed,DECOYS_for_TRAINING)

### entrenar glm segunda iteracion
glmFit_2<-train(Label~., data=train_DESCR_2[,c(2,4,5,7)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

### Pasar de nuevo los datos iniciales por el segundo glm

second_iteration<-predict(glmFit_2, TMP_OBJECT[,c(2,4,5,7)], type = 'prob')  #primera
TMP_OBJECT$p_Target_second<-as.numeric(second_iteration[,2])

SECOND_ITERATION_FDR<-psmFDR(TMP_OBJECT, "p_Target_second", "DECOY",concat_decoy=0, prot_id="Proteins")

SECOND_ITERATION_FDR_passed<-SECOND_ITERATION_FDR[which(SECOND_ITERATION_FDR$psmFDR <=0.01),]

############################# 3ª iteraciones
trainDESCR_3_tmp<-rbind(SECOND_ITERATION_FDR_passed[,c(1:19)], DECOYS_for_TRAINING)   #add decoys al objeto resultante de la iter anterior

tmp_train_third_iter_FDR<-psmFDR(trainDESCR_3_tmp, "hyperscore", "DECOY",concat_decoy=0, prot_id="Proteins")
tmp_train_third_iter_FDR_passed<-tmp_train_third_iter_FDR[which(tmp_train_third_iter_FDR$psmFDR<=0.01),]      #rerankea por FDR y quedate con lo que pase
tmp_train_third_iter_FDR_passed$psmFDR<-NULL

train_DESCR_3<-rbind(tmp_train_third_iter_FDR_passed,DECOYS_for_TRAINING) #añade de nuevo las decoys y crea objeto de entrenamiento

glmFit_3<-train(Label~., data=train_DESCR_3[,c(2,4,5,7)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

third_iteration<-predict(glmFit_3, TMP_OBJECT[,c(2,4,5,7)], type = 'prob')  #primera
TMP_OBJECT$p_Target_third<-as.numeric(third_iteration[,2])

THIRD_ITERATION_FDR<-psmFDR(TMP_OBJECT, "p_Target_third", "DECOY",concat_decoy=0, prot_id="Proteins")

THIRD_ITERATION_FDR_passed<-THIRD_ITERATION_FDR[which(THIRD_ITERATION_FDR$psmFDR <=0.01),]

#### 4ta

trainDESCR_4_tmp<-rbind(THIRD_ITERATION_FDR_passed[,c(1:19)], DECOYS_for_TRAINING)

tmp_train_fourth_iter_FDR<-psmFDR(trainDESCR_4_tmp, "hyperscore", "DECOY",concat_decoy=0, prot_id="Proteins")
tmp_train_fourth_iter_FDR_passed<-tmp_train_fourth_iter_FDR[which(tmp_train_fourth_iter_FDR$psmFDR<=0.01),]      #rerankea por FDR y quedate con lo que pase
tmp_train_fourth_iter_FDR_passed$psmFDR<-NULL

train_DESCR_4<-rbind(tmp_train_fourth_iter_FDR_passed,DECOYS_for_TRAINING) #añade de nuevo las decoys y crea objeto de entrenamiento

glmFit_4<-train(Label~., data=train_DESCR_4[,c(2,4,5,7)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

fourth_iteration<-predict(glmFit_4, TMP_OBJECT[,c(2,4,5,7)], type = 'prob')  #primera
TMP_OBJECT$p_Target_fourth<-as.numeric(fourth_iteration[,2])

FOURTH_ITERATION_FDR<-psmFDR(TMP_OBJECT, "p_Target_fourth", "DECOY",concat_decoy=0, prot_id="Proteins")
FOURTH_ITERATION_FDR_passed<-FOURTH_ITERATION_FDR[which(FOURTH_ITERATION_FDR$psmFDR <=0.01),]


###### 5ta
trainDESCR_5_tmp<-rbind(FOURTH_ITERATION_FDR_passed[,c(1:19)],DECOYS_for_TRAINING)

tmp_train_fifth_iter_FDR<-psmFDR(trainDESCR_5_tmp, "hyperscore", "DECOY",concat_decoy=0, prot_id="Proteins")
tmp_train_fifth_iter_FDR_passed<-tmp_train_fifth_iter_FDR[which(tmp_train_fifth_iter_FDR$psmFDR<=0.01),]      #rerankea por FDR y quedate con lo que pase
tmp_train_fifth_iter_FDR_passed$psmFDR<-NULL

train_DESCR_5<-rbind(tmp_train_fifth_iter_FDR_passed,DECOYS_for_TRAINING) #añade de nuevo las decoys y crea objeto de entrenamiento

glmFit_5<-train(Label~., data=train_DESCR_5[,c(2,4,5,7)], method="glm", family=binomial, trControl=ctrl, metric="ROC")

fifth_iteration<-predict(glmFit_5, TMP_OBJECT[,c(2,4,5,7)], type = 'prob')  #primera
TMP_OBJECT$p_Target_fifth<-as.numeric(fifth_iteration[,2])

FIFTH_ITERATION_FDR<-psmFDR(TMP_OBJECT, "p_Target_fifth", "DECOY",concat_decoy=0, prot_id="Proteins")
FIFTH_ITERATION_FDR_passed<-FIFTH_ITERATION_FDR[which(FIFTH_ITERATION_FDR$psmFDR <=0.01),]
