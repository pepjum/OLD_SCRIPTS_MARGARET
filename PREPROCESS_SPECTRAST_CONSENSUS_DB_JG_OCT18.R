#!/usr/bin/Rscript

args=(commandArgs(TRUE))

database<-args[1]   #sptxt de spectraST
spread_value<-args[2] #0   valor para máqunas modernas es 0
Dalton<-args[3]   #1
n_peaks<-args[4]   #filtering spectra by number of peaks
intensity<-args[5] #minimun intensity threshold
#type<-args[6]   #peaks or wavelets
wave<-args[6]  #dwt, dwtd, wp or peaks

#source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/margaret/data/pepe/scripts/functions_SE.R")

library(stringr)
library(stringi)
library(wavethresh)
library(compiler)
library(future.apply)
#library(feather)

options(stringsAsFactors = FALSE)
plan(multicore, workers=40)
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
list_all_peaks<-paste(db_spectra_filtered_filtered_by_intensity$Peaks)
db_spectra_filtered_filtered_by_intensity$Peaks<-NULL

cat("transform peaks info in dataframes to work in it","\n")

convert_df_2_c<-cmpfun(convert_df_2)
list_df_peaks_not_binned<-future_lapply(list_all_peaks, FUN = function(x) convert_df_2_c(x), future.seed=TRUE)


cat("bining","\n")
bining_c<-cmpfun(bining)
list_df_peaks_binned<-future_lapply(list_df_peaks_not_binned, FUN= function(x) bining_c(x,1), future.seed=TRUE)  #comprobar

if (spread_value !=0){
    cat("spread out","\n")
    list_df_peaks_binned_and_spread_out<-spread_out(list_df_peaks_binned,spread_value)
}else{
    list_df_peaks_binned_and_spread_out<-list_df_peaks_binned
}

cat("selecting range of spectra","\n")

select_range_c<-cmpfun(select_range)
list_selected<-future_lapply(list_df_peaks_binned_and_spread_out, FUN= function(x) select_range_c(x), future.seed=TRUE)

cat("normalizing peaks","\n")

normalization_c<-cmpfun(normalization)
list_df_normalized<-future_lapply(list_selected, FUN = function(x) normalization_c(x), future.seed=TRUE)
rm(list_selected)

cat("creating delimited vector","\n")

create_full_vector_c<-cmpfun(create_full_vector)
list_df_matrix<-future_lapply(list_df_normalized, FUN = function(x) create_full_vector_c(x), future.seed=TRUE)
rm(list_df_normalized)

if (wave == "peaks"){
    cat("transform dataframes into a list and rebuild  final database","\n")

    codification_c<-cmpfun(codification)
    list_linecoded<-future_lapply(list_df_matrix, FUN = function(x) codification_c(x), future.seed=TRUE)
    db_spectra_filtered_filtered_by_intensity$Peaks<-paste(list_linecoded)

    cat("saving peaks database","\n")
    new_name_p<-paste(substr(database, 1, nchar(database) - 5), "peaks", sep="")
    save(db_spectra_filtered_filtered_by_intensity, file=paste0(new_name_p,".Rdata"))
    numbers<-paste(db_spectra_filtered_filtered_by_intensity$PrecursorMZ)
    numbers_random<-sample(numbers)
    db_spectra_filtered_filtered_by_intensity$PrecursorMZ <- NULL
    db_spectra_filtered_filtered_by_intensity$PrecursorMZ <- paste(numbers_random)
    db_spectra_filtered_filtered_by_intensity$Name<-paste0(db_spectra_filtered_filtered_by_intensity$Name,"_DECOY")
    save(db_spectra_filtered_filtered_by_intensity, file=paste0(new_name_p, "_decoy.Rdata"))
}else{
    cat("transform to wavelets","\n")
    if(wave == "wp"){
    wavelets_c<-cmpfun(wavelets)
    wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)
    }else if(wave == "dwt"){
        wavelets_d_c<-cmpfun(wavelets_d)
        wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_d_c(x), future.seed=TRUE)
    }else if(wave == "dwtd"){
        wavelets_for_c<-cmpfun(wavelets_for)
        wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_for_c(x), future.seed=TRUE)

    }

    #normalization_w_c<-cmpfun(normalization_w)
    #wavelets_norm<-future_lapply(wavelets, FUN = function(x) normalization_w_c(x), future.seed=TRUE)
    cat("saving wavelets database","\n")
    new_name_w<-paste(substr(database, 1, nchar(database) - 5), "wavelets", sep="")
    db_spectra_filtered_filtered_by_intensity$Peaks<-NULL
    db_spectra_filtered_filtered_by_intensity$wavelet <- wavelets
    if(wave =="wp"){
        save(db_spectra_filtered_filtered_by_intensity, file=paste0(new_name_w,"_wp.Rdata"))
    }else if (wave =="dwt"){
        save(db_spectra_filtered_filtered_by_intensity, file=paste0(new_name_w,"_dwt.Rdata"))
    }else if (wave == "dwtd"){
        save(db_spectra_filtered_filtered_by_intensity, file=paste0(new_name_w,"_dwtd.Rdata"))
    }
    numbers<-paste(db_spectra_filtered_filtered_by_intensity$MW)
    numbers_random<-sample(numbers)
    db_spectra_filtered_filtered_by_intensity$MW <- NULL
    db_spectra_filtered_filtered_by_intensity$MW <- paste(numbers_random)
    db_spectra_filtered_filtered_by_intensity$Name<-paste0(db_spectra_filtered_filtered_by_intensity$Name,"_DECOY")
    if(wave =="wp"){
        save(db_spectra_filtered_filtered_by_intensity, file=paste0(new_name_w,"_decoy_wp.Rdata"))
    }else if (wave =="dwt"){
        save(db_spectra_filtered_filtered_by_intensity, file=paste0(new_name_w, "_decoy_dwt.Rdata"))
    }else if (wave == "dwtd"){
        save(db_spectra_filtered_filtered_by_intensity, file=paste0(new_name_w,"_decoy_dwtd.Rdata"))
    }

}
new<-Sys.time() - old
print(new)
cat("end of processing database","\n")
