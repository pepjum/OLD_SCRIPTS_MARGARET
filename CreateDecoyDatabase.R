args=(commandArgs(TRUE))
targetfastafile <- args[1]
decoyfastafile <- args[2]
decoyid <- args[3]

library(Biostrings)
source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

db_target <- readAAStringSet(targetfastafile, format="fasta",nrec=-1L, skip=0L, use.names=TRUE)
seq <- paste(db_target)
seq_decoy <- sapply(seq,  FUN=function(x) .strPseudoreverseDecoy(x))
names(seq_decoy) <- paste(decoyid,names(db_target),sep="_")
sp_decoy_fa <- AAStringSet(seq_decoy)
writeXStringSet(sp_decoy_fa, decoyfastafile)
