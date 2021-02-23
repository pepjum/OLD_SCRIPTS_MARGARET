args=(commandArgs(TRUE))


mgf<-args[1]      #mgf file
database<-args[2]   #database spectaST
rangeo<-args[3]     #for selecting spectra by mz difference of .... (0.1 by default)
spread_value<-args[4]    #spread out
Dalton<-args[5]       # Value of daltons for bining function purposes
output<-args[6]   # 1 for Rdata, 2 for txt, 3 for feather
best_number<-args[7]
#source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
source("/home/margaret/data/pepe/scripts/functions_SE.R")
#mgf_dir="."
cat("started at","\n")
old<-Sys.time()
print(old)

options(stringsAsFactors = FALSE)

library(stringr)
library(stringi)
library(plyr)
library(future.apply)
library(compiler)
#library(feather)
#setwd(paste0(mgf_dir,"/"))



enableJIT(3)
plan(multicore, workers=25)

options(future.globals.maxSize= +Inf)
cat("reading MGF file","\n")

ReadMGFFile_c<-cmpfun(ReadMGFFile)


mgf_file<-ReadMGFFile_c(mgf)


#READ SPECTRUMS AND KEEP NAMES AND M/Z CHARGE
mgf_to_df_c<-cmpfun(mgf_to_df4)
df_spectrums_mgf<-mgf_to_df(mgf_file)

rm(mgf_file)

cat("loading database","\n")
db_spectra<-get(load(database))   #_500-1500.Rdata

cat("selecting spectrums for MZ range","\n")
plan(multicore, workers=25)

options(future.globals.maxSize= +Inf)

selecting_spectrums_c<-cmpfun(selecting_spectrums)

df_spectrums_selected_for_each_spectrum<-ldply(future_apply(df_spectrums_mgf,1, function(x) selecting_spectrums_c(x,rangeo,db_spectra), future.seed=TRUE ),rbind)

list_all_peaks<-paste(df_spectrums_mgf$peaks)

cat("transform peaks info to a dataframe to work in it","\n")

convert_df_c<-cmpfun(convert_df)

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

cat("convert peaks dataframe into string to indexing in the info dataframe","\n")

codification_c<-cmpfun(codification)

list_linecoded<-future_lapply(list_df_matrix, FUN = function(x) codification_c(x), future.seed=TRUE)

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$Peaks_mgf<-paste(list_linecoded)

df_spectrums_mgf$idx<-seq(1:(nrow(df_spectrums_mgf)))

cat("merge","\n")
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)

cat("calculating dot product and dot bias","\n")

#save(df_results_tmp2,file="/home/margaret/data/pepe/por_si_fallo_dot_and_db.Rdata")

Dot_and_DB_comp<-cmpfun(Dot_and_DB)

res<-ldply(future_apply(df_results_tmp2,1, FUN=function(x) Dot_and_DB_comp(x), future.seed=TRUE),rbind)

res<-subset(res, select= -.id)

df_results_tmp2<-cbind(df_results_tmp2,res)


df_results_tmp2$Peaks_mgf<-NULL
df_results_tmp2$Peaks<- NULL
df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,-D)),]
df_results_tmp2<-df_results_tmp2[which(df_results_tmp2$D > 0),]


if (best_number !=0){
    cat("selecting only the best  ", best_number, "PSM","\n")
    selected_data<-list()
    for (i in unique(df_results_tmp2$idx)){
        tmp<-df_results_tmp2[which(df_results_tmp2$idx == i),]
        if(nrow(tmp)< best_number){
            tmp_selected<-tmp
        }else{
            tmp_selected<-tmp[1:best_number,]
        }
        tmp_selected <- na.omit(tmp_selected)
        selected_data[[i]]<-tmp_selected
    }
    df_results_tmp2<-ldply(selected_data,rbind)
}

df_results_tmp2[which(is.nan(df_results_tmp2$DB)), "DB"]<-0


cat("Calculating delta","\n")

Delta_c<-cmpfun(Delta)

df_results_tmp3<-Delta_c(df_results_tmp2)
rm(df_results_tmp2)

cat("Fscore","\n")

F_score_com<-cmpfun(F_score)

res_2<-future_apply(df_results_tmp3,1, FUN = function(x) F_score(x), future.seed=TRUE)
res_2<-as.data.frame(res_2)
names(res_2)<-"Fscore"

df_results_tmp2<-cbind(df_results_tmp3,res_2)
new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")

if(output == 1){
    save(df_results_tmp2, file = paste0(new_name,".SpectraST.Rdata"))
}else if(output == 2){
    write.table(df_results_tmp2, file= paste0(new_name,".txt"), row.names=F, col.names=F, quote = FALSE, sep="\t")
}

new<-Sys.time() - old
print(new)
cat("end of search")
