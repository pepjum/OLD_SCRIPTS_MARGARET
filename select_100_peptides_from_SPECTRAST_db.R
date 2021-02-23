
source("/home/margaret/data/pepe/scripts/functions_SE.R")

library(stringr)
library(stringi)
library(wavethresh)
library(compiler)
library(future.apply)
#library(feather)

options(stringsAsFactors = FALSE)
plan(multicore, workers=25)
enableJIT(3)
options(future.globals.maxSize= +Inf)

old<-Sys.time()
print(old)

#####CALL TO DB
cat("loading db","\n")
cat("Patience.... ","\n")

read_Spectra_C<-cmpfun(readSpectraSTdb)
old<-Sys.time()
db_spectra<-read_Spectra_C(paste(database))

cat("filtering spectra by number of peaks","\n")

filteringSpectraByNumberOfPeaks_c<-cmpfun(filteringSpectraByNumberOfPeaks)
db_spectra_filtered<-filteringSpectraByNumberOfPeaks_c(db_spectra,n_peaks)

cat("filtering spectra by intensity","\n")

filteringSpectraByIntensity_c<-cmpfun(filteringSpectraByIntensity)
db_spectra_filtered_filtered_by_intensity<-filteringSpectraByIntensity_c(db_spectra_filtered, intensity)


output_list<-c()
for(i in 1:nrow(db_spectra_filtered_filtered_by_intensity)){
    spec<-db_spectra_filtered_filtered_by_intensity[i,]
    a<-"BEGIN IONS"
    output_list<-c(output_list,a)
    Title<-paste0("TITLE=",spec$scan)
    output_list<-c(output_list,Title)
    Pepmass<-paste0("PEPMASS=", spec$mass/2)
    output_list<-c(output_list,Pepmass)
    Charge<-paste0("CHARGE=",spec$charge)
    output_list<-c(output_list,Charge)
    Seqs<-paste0("SEQ=",spec$peptides)
    output_list<-c(output_list,Seqs)
    peaks<-as.data.frame(list_dataframes[i])
    for (k in 1:nrow(peaks)){
        lines<-peaks[k,]
        line<-paste(lines$m,lines$z, sep="\t")
        output_list<-c(output_list,line)
    }
    b<-"END IONS\n"
    output_list<-c(output_list,b)
}

writeLines(unlist(lapply(output_list, paste, collapse="\n")), "/home/margaret/data/pepe/100_peptides_from_SPECTRAST_DB.mgf")
