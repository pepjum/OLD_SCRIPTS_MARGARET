source("/home/margaret/data/pepe/scripts/functions_SE.R")
#mgf_dir="."
cat("started at","\n")
old<-Sys.time()
print(old)

options(stringsAsFactors = FALSE)

mgf<-"/home/margaret/data/pepe/100_peptides_from_SPECTRAST_DB.mgf"


library(stringr)
library(stringi)
library(plyr)
library(future.apply)
library(compiler)
library(wavethresh)
#setwd(paste0(mgf_dir,"/"))

enableJIT(3)
plan(multicore, workers=25)

options(future.globals.maxSize= +Inf)
cat("reading MGF file","\n")

ReadMGFFile_c<-cmpfun(ReadMGFFileMASSSimulator)


mgf_file<-ReadMGFFile_c(mgf)


#READ SPECTRUMS AND KEEP NAMES AND M/Z CHARGE
mgf_to_df_c<-cmpfun(mgf_to_df4)
#df_spectrums_mgf<-mgf_to_df(mgf_file)
df_spectrums_mgf<-mgf_to_df_c(mgf_file)


list_all_peaks<-paste(df_spectrums_mgf$peaks)

cat("transform peaks info to a dataframe to work in it","\n")

convert_df_c<-cmpfun(convert_df3)

list_df_peaks_not_binned<-future_lapply(list_all_peaks, FUN = function(x) convert_df_c(x), future.seed=TRUE)


cat("bining","\n")

bining_c<-cmpfun(bining)
list_df_peaks_binned<-future_lapply(list_df_peaks_not_binned, FUN= function(x) bining_c(x,1), future.seed=TRUE)


if (spread_value !=0){
    cat("spread_out","\n")
    list_df_peaks_binned_and_spread_out<-spread_out(list_df_peaks_binned,spread_value)
}else{
    list_df_peaks_binned_and_spread_out<-list_df_peaks_binned
}

cat("selecting range of peak list","\n")

select_range_c<-cmpfun(select_range)
list_selected<-future_lapply(list_df_peaks_binned_and_spread_out, FUN= function(x) select_range_c(x), future.seed=TRUE)

cat("normalization","\n")
normalization_c<-cmpfun(normalization)

list_df_normalized<-future_lapply(list_selected, FUN = function(x) normalization_c(x), future.seed=TRUE)

cat("creating a vector adding zeros to positions where we have not mz values","\n")

create_full_vector_c<-cmpfun(create_full_vector)

list_df_matrix<-future_lapply(list_df_normalized, FUN = function(x) create_full_vector_c(x), future.seed=TRUE)

mu<-c()
sd<-c()
for(i in 1:length(list_df_matrix)){
    spec<-list_df_matrix[[i]]
    sd_spec<-sd(spec$normalized)
    sd<-c(sd,sd_spec)
    mu_spec<-mean(spec$normalized)
    mu<-c(mu,mu_spec)
}

media_sd<-mean(sd)  #0.02697347
media_mu<-mean(mu) #0.003563984
