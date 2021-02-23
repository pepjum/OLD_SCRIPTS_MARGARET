### PREPARAR datos en R



library(seqinr)

library(protr) # protr::AAindex (544 properties) podría intentar usarse con la función de Interpol)
library(Peptides)

data<-get(load("~/data/pepe/18_CLASIFICADOR_DP_OCT19/MATRICES/ALL_EXPERIMENT_WITHOUT_AMBIGUOUS_TyD_PREPARED.Rdata"))

data$PeptideMass<-NULL


data$PeptideMass<-mw(paste(data$PeptideSeq), monoisotopic = TRUE)


write.table(data,file="ALL_EXPERIMENT_WITHOUT_AMBIGUOUS_TyD_PREPARED.txt", col.names=F, row.names=F, sep="\t", quote=F)
