setwd("/home/bioinformatica/datos/03_Analysis/eguruce/16_aaindex_Mar16")

# save.image("classifier_abr16.RData")
# load("classifier_abr16.RData")
library(beeswarm)
library(seqinr)
data(aaindex)
library(e1071)
library(ROCR)
# library(pROC)
# library(rJava)

# load("processedPepInput.rda")
# load("processedPepNObsInput.rda")
# load("all_prop.rda")
# load("sample_500.rda")
# load("sample_500_Control.rda")


##################################################################################################################################################
# 1. Selección de las propiedades con FDR<0.05 en todos los test del bootstrapping-rarefaction (MissingObs vs NoMissing)
# 1b. Selección de las propiedades con FDR<0.05 en todos los test del bootstrapping-rarefaction (MissingNObs vs NoMissing): No hay ninguna
# 2. Comprobar que no tienen p<0.05 en el control - (NoMissing vs NoMissing)
# 3. Escoger la propiedad representante de cada grupo de propiedades correladas (utilizo la info de correlación de la DB aaindex) que tenga < FDR
# 4. Selección de péptidos para el training y el test y propiedades óptimas para el clasificador (caret)
# 5. Generar el clasificador y validarlo con los péptidos seleccionados para el test en otro script
#--------- Esta selección de péptidos para el clasificador no vale. Ir a "aaindex_gpmdb_EG_may16.R"

aaindex_annot <- data.frame(ID=unlist(lapply(aaindex, FUN=function(x) unlist(x)["H"])), Property=unlist(lapply(aaindex, FUN=function(x) unlist(x)["D"])), Reference=unlist(lapply(aaindex, FUN=function(x) unlist(x)["R"])), Authors=unlist(lapply(aaindex, FUN=function(x) unlist(x)["A"])), Description=unlist(lapply(aaindex, FUN=function(x) unlist(x)["T"])), Journal=unlist(lapply(aaindex, FUN=function(x) unlist(x)["J"])), Cor=unlist(lapply(aaindex, FUN=function(x) unlist(x)["C"])))
rownames(aaindex_annot) <- aaindex_annot[,1]

#--------------------------------------------------------------------------------------------------------------------------------------------------
# PROPERTIES per Peptide
# "For each peptide and a given property, the constituent amino acid numerical values were averaged to produce a single value. Missing values were ignored. The average (rather than median or sum) was chosen because it is sensitive to outliers and normalizes for peptide length. It was assumed that the average physicochemical property across each peptide was sufficient to capture relevant information about peptide response." Fusaro et al. 2009

# 500 iteraciones con el 80 % de los péptidos missing con observaciones vs igual número de péptidos no-missing (rownames->prop , colnames=iteration) - 2490 genes
for(i in 1:10)
{
	load(paste("cluster/rarefaction/",i*50-49,"DB1_rarefaction.rda",sep=""))
	tmp <- rarefaction_pvalue
	for(j in 2:11)
	{
		load(paste("cluster/rarefaction/",i*50-49,"DB",j,"_rarefaction.rda",sep=""))
		tmp <- rbind(tmp, rarefaction_pvalue)
	}
	if(i==1) rarefaction_pvalue_cluster <- tmp else rarefaction_pvalue_cluster <- cbind(rarefaction_pvalue_cluster,tmp)
}
rarefaction_fdr_cluster <- apply(rarefaction_pvalue_cluster, 2, FUN=function(x) p.adjust(x, method="fdr"))
sum(apply(rarefaction_fdr_cluster, 1, FUN=function(x) sum(x<0.05))==500) # 293 con FDR<0.05 en todos los muestreos - 261 con FDR<0.01

prop_sig <- names(which(apply(rarefaction_fdr_cluster, 1, FUN=function(x) sum(x<0.05))==500))
aaindex_annot_sig <- aaindex_annot[prop_sig[-c(286:293)],][which(aaindex_annot[prop_sig[-c(286:293)],"Cor"]==""),]
# Bucle para quedarse con la propiedad representativa de un conjunto de propiedades correlacionadas. Selecciono la que tiene < FDR
prop_sig_cor <- setdiff(prop_sig[-c(286:293)],rownames(aaindex_annot_sig))
prop_sig_cor_tmp <- data.frame()
prop_grp_all <- vector()
for(i in 1:length(prop_sig_cor))
{
	if(i==1)
	{
		prop_grp <- c(prop_sig_cor[i],intersect(unique(unlist(strsplit(gsub("[0-9][.][0-9][0-9][0-9]","",paste(aaindex_annot[prop_sig_cor[i],"Cor"]))," "))),prop_sig_cor))
		prop_sig_cor_tmp <- aaindex_annot[names(which.min(sapply(prop_grp, FUN=function(x) mean(rarefaction_fdr_cluster[x,])))),]
		prop_grp_all <- c(prop_grp_all,prop_grp)
		print(i)
	}else{
		if(sum(prop_sig_cor[i]==prop_grp_all)==0)
		{
			prop_grp <- c(prop_sig_cor[i],intersect(unique(unlist(strsplit(gsub("[0-9][.][0-9][0-9][0-9]","",paste(aaindex_annot[prop_sig_cor[i],"Cor"]))," "))),prop_sig_cor))
			prop_sig_cor_tmp <- rbind(prop_sig_cor_tmp, aaindex_annot[names(which.min(sapply(prop_grp, FUN=function(x) mean(rarefaction_fdr_cluster[x,])))),])
			prop_grp_all <- c(prop_grp_all,prop_grp)
			print(i)
		}
	}
}
aaindex_annot_sig_all <- rbind(aaindex_annot_sig,unique(prop_sig_cor_tmp), cbind(ID=prop_sig[c(286:293)], Property=prop_sig[c(286:293)], Reference=rep("seqinr",length(prop_sig[c(286:293)])),Authors="",Description="",Journal="",Cor=""))
aaindex_annot_sig_all$Average_fdr <- apply(rarefaction_fdr_cluster[paste(aaindex_annot_sig_all[,1]),], 1, mean)

# 500 iteraciones con el 80 % de los péptidos missing sin observaciones vs igual número de péptidos no-missing con observaciones (rownames->prop , colnames=iteration) - 201138 genes
for(i in 1:10)
{
  load(paste("cluster/rarefaction/",i*50-49,"DB1_nobs_rarefaction.rda",sep=""))
  tmp <- rarefaction_nobs_pvalue
  for(j in 2:11)
  {
    load(paste("cluster/rarefaction/",i*50-49,"DB",j,"_nobs_rarefaction.rda",sep=""))
    tmp <- rbind(tmp, rarefaction_nobs_pvalue)
  }
  if(i==1) rarefaction_nobs_pvalue_cluster <- tmp else rarefaction_nobs_pvalue_cluster <- cbind(rarefaction_nobs_pvalue_cluster,tmp)
}
rarefaction_nobs_fdr_cluster <- apply(rarefaction_nobs_pvalue_cluster, 2, FUN=function(x) p.adjust(x, method="fdr"))
sum(apply(rarefaction_nobs_pvalue_cluster, 1, FUN=function(x) sum(x<0.05))==500) # 293 con FDR<0.05 en todos los muestreos - 261 con FDR<0.01
aaindex_annot[names(which(apply(rarefaction_nobs_pvalue_cluster, 1, FUN=function(x) sum(x<0.05))==10)), ]


# Control -: 500 iteraciones con mismo número de péptidos que en el análisis de las missing pero comparando no-missing vs no-missing (rownames->prop , colnames=iteration)
for(i in 1:10)
{
	load(paste("cluster/rarefaction/",i*50-49,"DB1_rarefaction_C.rda",sep=""))
	tmp <- rarefaction_pvalue_controls
	for(j in 2:11)
	{
		load(paste("cluster/rarefaction/",i*50-49,"DB",j,"_rarefaction_C.rda",sep=""))
		tmp <- rbind(tmp, rarefaction_pvalue_controls)
	}
	if(i==1) rarefaction_pvalue_control_c <- tmp else rarefaction_pvalue_control_c <- cbind(rarefaction_pvalue_control_c,tmp)
}
rarefaction_fdr_control_c <- apply(rarefaction_pvalue_control_c, 2, FUN=function(x) p.adjust(x, method="fdr"))
sum(apply(rarefaction_pvalue_control_c, 1, FUN=function(x) sum(x<0.05))==500) # 0 con p<0.05 en todos los muestreos


pdf(file="hist_p_fdr_may16.pdf")
{
par(mfrow=c(3,2))
hist(rarefaction_pvalue_cluster, main="Missing vs Non-Missing (p-value)")
hist(rarefaction_fdr_cluster, main="Missing vs Non-Missing (FDR)")
hist(rarefaction_nobs_pvalue_cluster, main="Missing_nobs vs Non-Missing (p-value)")
hist(rarefaction_nobs_fdr_cluster, main="Missing_nobs vs Non-Missing (FDR)")
hist(rarefaction_pvalue_control_c, main="Non-Missing vs Non-Missing (p-value)")
hist(rarefaction_fdr_control_c, main="Non-Missing vs Non-Missing (FDR)")
}
dev.off()

#----------------------------------------------------------------------------------------------------------------------------------------------------
# CLASSIFIER

# number of observations of the peptide sequence with a parent ion charge of +1 (z.1), +2, +3 or +4 and minimum observed E-value (E >= 0.01 means peptide not observed EC=1) http://wiki.thegpm.org/wiki/GPMDB_evidence_codes
pdf(file="hist_gpmdb_may16.pdf")
{
par(mfrow=c(2,2))
hist(log10(pepObs_gpmdb_proteotypicXmis_f2$Obs), main="GPMDB log10(Obs)")
hist(-log10(pepObs_gpmdb_proteotypicXmis_f2$E), main="GPMDB -log10(E)")
hist(log10(pepObs_gpmdb_proteotypicXmis_f2$Obs)[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1], main="GPMDB missing log10(Obs)")
hist(log10(pepObs_gpmdb_proteotypicXmis_f2$Obs)[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0], main="GPMDB non-missing log10(Obs)")
}
dev.off()
pdf(file="hist2_gpmdb_may16.pdf")
{
hist(log10(pepObs_gpmdb_proteotypicXmis_f2$Obs)[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0],100, main="GPMDB missing log10(Obs)", xlim=c(0,8))
}
dev.off()

pepObs_gpmdb_proteotypicXmis_f2$Obs <- rowSums(pepObs_gpmdb_proteotypicXmis_f2[,2:5])
GPMDBe_MissVsNoMiss_t <- t.test(pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1,"E"],pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0,"E"], na.rm=TRUE)
GPMDBobs_MissVsNoMiss_t <- t.test(pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==1,"Obs"],pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2$Missing.min==0,"Obs"], na.rm=TRUE)
pdf("Missing_E.GPMDB_Boxplot_May16.pdf", width = 10, height = 10)
{

	d <- data.frame(expression=-log10(pepObs_gpmdb_proteotypicXmis_f2[,"E"]),lab=factor(pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"], label=c("Non-Missing","Missing")))

	boxplot(expression ~ lab, data = d, outline = FALSE, xlab = " ", ylab = "-log10(E-value)", main = paste("Missing vs Non-Missing (p-value=",round(GPMDBe_MissVsNoMiss_t$p.value, 70),")", sep=""))

	# beeswarm(expression ~ lab, data = d, pch = 21, col = c("blue","orange"), bg = "#00000050", add = TRUE, method="swarm", corral = "random")
}
dev.off()
pdf("Missing_Obs.GPMDB_Boxplot_May16.pdf", width = 10, height = 10)
{

	d <- data.frame(expression=log10(pepObs_gpmdb_proteotypicXmis_f2[,"Obs"]),lab=factor(pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"], label=c("Non-Missing","Missing")))

	boxplot(expression ~ lab, data = d, outline = FALSE, xlab = " ", ylab = "log10(Observations)", main = paste("Missing vs Non-Missing (p-value=",round(GPMDBobs_MissVsNoMiss_t$p.value, 70),")", sep=""))

	# beeswarm(expression ~ lab, data = d, pch = 21, col = c("blue","orange"), bg = "#00000050", add = TRUE, method="swarm", corral = "random")
}
dev.off()

# Peptide selection
pdf(file = "ExObs_May16.pdf", width = 8, height = 8)
plot(log10(pepObs_gpmdb_proteotypicXmis_f2$Obs), -log10(pepObs_gpmdb_proteotypicXmis_f2$E), xlab="log10(Obs)", ylab="-log10(E-value)", pch=".")
dev.off()

pepObs_gpmdb_proteotypicXmis_f2 <- pepObs_gpmdb_proteotypicXmis_f2[order(pepObs_gpmdb_proteotypicXmis_f2$E, decreasing=FALSE),]
rownames(pepObs_gpmdb_proteotypicXmis_f2) <- pepObs_gpmdb_proteotypicXmis_f2[,1]
trainE_pep <- rbind(pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"]==0, ][sample(2000,1500),], pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"]==1, ][sample(2000,1500),])
testE_pep <- rbind(pepObs_gpmdb_proteotypicXmis_f2[setdiff(paste(pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"]==0, ][1:2000,1]),paste(trainE_pep[1:1500,1])),], pepObs_gpmdb_proteotypicXmis_f2[setdiff(paste(pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"]==1, ][1:2000,1]),paste(trainE_pep[1501:3000,1])),])

pepObs_gpmdb_proteotypicXmis_f2 <- pepObs_gpmdb_proteotypicXmis_f2[order(pepObs_gpmdb_proteotypicXmis_f2$Obs, decreasing=TRUE),]
trainObs_pep <- rbind(pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"]==0, ][sample(2000,1500),], pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"]==1, ][sample(2000,1500),])
testObs_pep <- rbind(pepObs_gpmdb_proteotypicXmis_f2[setdiff(paste(pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"]==0, ][1:2000,1]),paste(trainObs_pep[1:1500,1])),], pepObs_gpmdb_proteotypicXmis_f2[setdiff(paste(pepObs_gpmdb_proteotypicXmis_f2[pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"]==1, ][1:2000,1]),paste(trainObs_pep[1501:3000,1])),])

#---------------------------------
## NAIVE BAYES
## El paquete e1071 no funciona bien para una variable cuantitativa así que me paso a klaR

trainE_pep
trainObs_pep
aaindex_annot_sig_all

matClassE <- data.frame(apply(data.frame(trainE_pep[, c("sequence", "Missing.min")], t(all_prop[rownames(aaindex_annot_sig_all),rownames(trainE_pep)])), 2, FUN=function(x) unlist(x)))
matClassO <- data.frame(apply(data.frame(trainObs_pep[, c("sequence", "Missing.min")], t(all_prop[rownames(aaindex_annot_sig_all),rownames(trainObs_pep)])), 2, FUN=function(x) unlist(x)))

a <- data.frame(Length=(as.numeric(paste(all_prop["Length",paste(pepObs_gpmdb_proteotypicXmis_f2[,1])]))>12)*1)
b <- factor(pepObs_gpmdb_proteotypicXmis_f2[,"Missing.min"], levels=c("1","0"))
classifier <- naiveBayes(as.numeric(paste(all_prop["Length",paste(pepObs_gpmdb_proteotypicXmis_f2[,1])])), b)
probpred <- predict(classifier, as.numeric(paste(all_prop["Length",paste(pepObs_gpmdb_proteotypicXmis_f2[,1])])), type = "raw")
predictMissing <- probpred[,1]
classifer_k <- NaiveBayes(as.numeric(paste(all_prop["Length",paste(pepObs_gpmdb_proteotypicXmis_f2[,1])])), b)
probpred_k <- predict(classifer_k)

classifier2 <- naiveBayes(a, b)
probpred2 <- predict(classifier,a, type = "raw")
predictMissing2 <- probpred[,1]

classifier_sim <- naiveBayes(data.frame(L=c(rep(1,15), c(3:15,3:15)),L2=c(rep(1,15), c(3:15,3:15))), factor(c(rep(1,15),rep(0,26)), levels=c("1","0")))
classifier_sim <- naiveBayes(c(rep(1,15), c(3:15,3:15)), factor(c(rep(1,15),rep(0,26)), levels=c("1","0")))
probpred_sim <- predict(classifier_sim,c(rep(1,15), c(3:15,3:15)), type = "raw")
probpred_sim <- predict(classifier_sim,data.frame(L=c(rep(1,15), c(3:15,3:15)),L2=c(c(rep(1,15), c(3:15,3:15)))), type = "raw")
predictMissing_sim <- probpred_sim[,1]

pred <- prediction(predictMissing, b)
perf <- performance(pred, "tpr", "fpr")
plot(perf)
performance(pred, "auc")@y.values

# Leer funciones del final del script para gráfica con ggplot2
roc_missing <- rocdata(b, predictMissing)
p <- rocplot.single(roc_missing, title = "", TRUE)
predPP <- rbind(data.frame("Pep" = paste(matClassE[matClassE[,2] == 1,1]), "Set" = "Missing", "Posterior Probability" = predictMissing[matClassE[,2] == 1]), data.frame("Pep" = paste(matClassE[matClassE[,2] == 0,1]), "Set" = "Known", "Posterior Probability" = predictMissing[matClassE[,2] == 0]))

pdf(file = "NaiveBayesProp1.2_May16.pdf", width = 8, height = 8)
ggplot(predPP, aes(factor(Set), Posterior.Probability)) + geom_boxplot(fill = "grey80", colour = "#3366FF")
dev.off()
pdf(file = "NaiveBayesProp1ROC_May16.pdf", width = 8, height = 8)
p
dev.off()
score_E <- data.frame("Seq" = paste(matClassE[,1]), "ProbLength" = predictMissing)

library(klaR)
locpvs
NaiveBayes
#---------------------------------
## NEURONAL network
# Interpol
library(randomForest)
library(ROCR)
rf <- randomForest(as.factor(classes),data = hydrophobicity_norm_l10) # build forest
pred <- prediction(rf$votes[,2], classes)             #prediction object
perf <- performance(pred, "auc")

# protR
rf.fit = randomForest(x.tr, y.tr, cv.fold = 5)
print(rf.fit)
rf.pred = predict(rf.fit, newdata = x.te, type ='prob')[, 1]
plot.roc(y.te, rf.pred, col ='#0080ff', grid = TRUE, print.auc = TRUE)




#####################################################################################
##
## FUNCIONES ROC
##
#####################################################################################

"rocdata" <- function(grp, pred){
  # Produces x and y co-ordinates for ROC curve plot
  # Arguments: grp - labels classifying subject status
  #            pred - values of each observation
  # Output: List with 2 components:
  #         roc = data.frame with x and y co-ordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval

  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  }

  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }

  cut <- unique(pred)
  tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]

  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2

  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)

  se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
  z <- (auc - 0.5)/se.auc.null
  p <- 2*pnorm(-abs(z))

  stats <- data.frame (auc = auc,
                       p.value = p,
                       ci.upper = ci.upper,
                       ci.lower = ci.lower
                       )

  return (list(roc = roc, stats = stats))
}

"rocplot.single" <- function(rocdata, title = "ROC Plot", p.value = FALSE){
  require(ggplot2)
  plotdata <- rocdata

  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
  }

  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
      geom_line(aes(colour = "")) +
      geom_abline (intercept = 0, slope = 1) +
      theme_bw() +
      scale_x_continuous("False Positive Rate (1-Specificity)") +
      scale_y_continuous("True Positive Rate (Sensitivity)") +
      scale_colour_manual(labels = annotation, values = "#000000") +
      theme(title = element_text(title, face="bold", size=14),
            #plot.title = theme_text(),
           axis.title.x = element_text(face="bold", size=12),
           axis.title.y = element_text(face="bold", size=12, angle=90),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.justification=c(1,0),
           legend.position=c(1,0),
           legend.title=element_blank(),
           legend.key = element_blank()
           )
  return(p)
}

"rocplot.multiple" <- function(test.data.list, groupName = "grp", predName = "res", title = "ROC Plot", p.value = TRUE){
  require(plyr)
  require(ggplot2)
  plotdata <- llply(test.data.list, function(x) with(x, rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
  plotdata <- list(roc = ldply(plotdata, function(x) x$roc),
                   stats = ldply(plotdata, function(x) x$stats)
                   )

  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
  }

  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
       geom_line(aes(colour = .id)) +
       geom_abline (intercept = 0, slope = 1) +
       theme_bw() +
       scale_x_continuous("False Positive Rate (1-Specificity)") +
       scale_y_continuous("True Positive Rate (Sensitivity)") +
       scale_colour_brewer(palette="Set1", breaks = names(test.data.list), labels = paste(names(test.data.list), ": ", annotation, sep = "")) +
       theme(title = title,
            plot.title = theme_text(face="bold", size=14),
            axis.title.x = theme_text(face="bold", size=12),
            axis.title.y = theme_text(face="bold", size=12, angle=90),
            panel.grid.major = theme_blank(),
            panel.grid.minor = theme_blank(),
            legend.justification=c(1,0),
            legend.position=c(1,0),
            legend.title=theme_blank(),
            legend.key = theme_blank()
            )
  return(p)
}
