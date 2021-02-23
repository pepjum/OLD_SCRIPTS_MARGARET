args=(commandArgs(TRUE))


mgf<<-args[1]      #mgf file
database<<-args[2]   #database spectaST
rango<<-args[3]     #for selecting spectra by mz difference of .... (0.1 by default)
spread_value<<-args[4]    #spread out
Dalton<<-args[5]       # Value of daltons for bining function purposes
output<-args[6]   # 1 for Rdata, 2 for txt, 3 for feather


mgf<<-"/home/margaret/data/pepe/PXD001381/21877.mgf"
database<<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.sptxt_range500-1500.Rdata"
rango<<-0.5     #for selecting spectra by mz difference of .... (0.1 by default)
spread_value<<-0    #spread out
Dalton<<-1       # Value of daltons for bining function purposes



source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")
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
library(feather)
#library(foreach)
#library(doParallel)

#setwd(paste0(mgf_dir,"/"))

enableJIT(3)
cat("reading MGF file","\n")

ReadMGFFile_c<-cmpfun(ReadMGFFile)
convert_df_c<-cmpfun(convert_df)
selecting_spectrums_c<-cmpfun(selecting_spectrums)
bining_c<-cmpfun(bining)
select_range_c<-cmpfun(select_rango)
normalization_c<-cmpfun(normalization)
create_full_vector_c<-cmpfun(create_full_vector)
codification_c<-cmpfun(codification)
Dot_and_DB_comp<-cmpfun(Dot_and_DB)
Delta_c<-cmpfun(Delta)
F_score_com<-cmpfun(F_score)





#cl<-makeCluster(25)
#registerDoParallel(cl)

mgf_file<-ReadMGFFile_c(mgf)

name_spectrums<-c()
mass_values<-c()
peaks_list<-c()

#READ SPECTRUMS AND KEEP NAMES AND M/Z CHARGE

for (i in 1:length(mgf_file)){
    name_file<-mgf_file[[i]][2]
    name_spectrums<-c(name_spectrums,name_file)
    mass<-mgf_file[[i]][4]
    mass_values<-c(mass_values,mass)
    peaks<-mgf_file[[i]][6:length(mgf_file[[i]])-1]
    peaks_cod<-paste(peaks,collapse="###")
    peaks_list<-c(peaks_list, peaks_cod)

}
mass_values<-sapply(strsplit(paste(mass_values),"="),`[`,2)
mass_values<-sapply(strsplit(paste(mass_values)," "),`[`,1)

name_files<-sapply(strsplit(paste(name_spectrums),"="),`[`,2)
name_spectrums<-sapply(strsplit(paste(name_files)," "),`[`,1)

#creating df from mgf file

df_spectrums_mgf<-data.frame("spectrum"=name_spectrums,"mz_value"=mass_values, "peaks"=paste(peaks_list))
df_spectrums_mgf$idx<-seq(1:nrow(df_spectrums_mgf))

rm(mgf_file)

cat("loading database","\n")
db_spectra<-get(load(database))   #_500-1500.Rdata

plan(multicore, workers=25)


#list_df_for_each_spec<-list()
cat("searching....","\n")
#for (i in 1:nrow(df_spectrums_mgf)){
### Define as global ,range,spread_value,Dalton,db_spectra
search_mgf<-function(dataframe_mgf){
    # spec<-as.data.frame(dataframe_mgf)
    dataframe_mgf <- as.data.frame(dataframe_mgf)
    options(future.globals.maxSize= +Inf)
    print(dataframe_mgf[1,1])
    #cat(paste(i, ifelse(i %% 100 == 0,"\n"," ")))
    #cat("+++++++++++++++++++++++++++++")
    # df_spectrums_selected_for_each_spectrum<-ldply(apply(dataframe_mgf,1, function(x) selecting_spectrums_c(x,rango,db_spectra)),rbind)
    print(class(dataframe_mgf))
    # selecting_spectrumsManuel<-function(df_spectrums_mgf,rango,db_spectra){
    #         PrecursorMGF<-as.numeric(paste(df_spectrums_mgf$mz_value))
    #         i<-as.numeric(paste(df_spectrums_mgf$idx))
    #         last<-selectSpectraByMZrango(db_spectra,PrecursorMGF,rango, i)
    # }
    #
    # df_spectrums_selected_for_each_spectrum <- selecting_spectrumsManuel(dataframe_mgf, rango, db_spectra)
    PrecursorMGF<-as.numeric(paste(dataframe_mgf$mz_value))
    i<-as.numeric(paste(df_spectrums_mgf$idx))
    df_spectrums_selected_for_each_spectrum <-selectSpectraByMZrango(db_spectra,PrecursorMGF,rango, i)
    print(class(df_spectrums_selected_for_each_spectrum))
    print(paste0('checking...... ', df_spectrums_selected_for_each_spectrum[1,1]))
    spec_rows <- nrow(df_spectrums_selected_for_each_spectrum)
    print(spec_rows)
    ifelse(spec_rows !=0,matchingSpec <- T,matchingSpec <- F)
    print(matchingSpec)
    if (matchingSpec){
        print('inside')
        list_all_peaks<-paste(dataframe_mgf$peaks)
        #cat("transform peaks info to a dataframe to work in it","\n")
        list_df_peaks_not_binned<-future_lapply(list_all_peaks, FUN = function(x) convert_df_c(x), future.seed=TRUE)
        #cat("bining","\n")
        list_df_peaks_binned<-future_lapply(list_df_peaks_not_binned, FUN= function(x) bining_c(x,Dalton), future.seed=TRUE)

        if (spread_value !=0){
            cat("spread_out","\n")
            list_df_peaks_binned_and_spread_out<-spread_out(list_df_peaks_binned,spread_value)
        }else{
            list_df_peaks_binned_and_spread_out<-list_df_peaks_binned
        }

        #cat("selecting range of peak list","\n")

        list_selected<-future_lapply(list_df_peaks_binned_and_spread_out, FUN= function(x) select_range_c(x), future.seed=TRUE)
        #cat("normalization","\n")

        list_df_normalized<-future_lapply(list_selected, FUN = function(x) normalization_c(x), future.seed=TRUE)

        #cat("creating a vector adding zeros to positions where we have not mz values","\n")

        list_df_matrix<-future_lapply(list_df_normalized, FUN = function(x) create_full_vector_c(x), future.seed=TRUE)

        #cat("convert peaks dataframe into string to indexing in the info dataframe","\n")

        list_linecoded<-future_lapply(list_df_matrix, FUN = function(x) codification_c(x), future.seed=TRUE)

        dataframe_mgf$peaks<-NULL
        dataframe_mgf$Peaks_mgf<-paste(list_linecoded)

        #cat("merge","\n")
        df_results_tmp<-merge(dataframe_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )
        df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
        rm(df_results_tmp)

        if(nrow(df_results_tmp2 !=0)){
            #cat("calculating dot product and dot bias","\n")



            res<-ldply(future_apply(df_results_tmp2,1, FUN=function(x) Dot_and_DB_comp(x), future.seed=TRUE), rbind)
            res<-subset(res, select= -.id)

            df_results_tmp2<-cbind(df_results_tmp2,res)
            df_results_tmp2$Peaks_mgf<-NULL
            df_results_tmp2$Peaks<- NULL
            df_results_tmp2<-df_results_tmp2[which(df_results_tmp2$D > 0),]
            df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,-D)),]

            #cat("calculating DeltaD","\n")
            df_results_tmp2$Delta<-NA
            df_results_tmp2[which(is.nan(df_results_tmp2$DB)), "DB"]<-0

            df_results_tmp3<-Delta_c(df_results_tmp2)
            rm(df_results_tmp2)

            #cat("calculating Fscore","\n")

            res_2<-future_apply(df_results_tmp3,1, FUN = function(x) F_score(x), future.seed=TRUE)
            res_2<-as.data.frame(res_2)
            names(res_2)<-"Fscore"
            df_results_tmp2<-cbind(df_results_tmp3,res_2)
            return(df_results_tmp2)
            #list_df_for_each_spec[[i]]<-df_results_tmp2
        }
    }else{return(NA)}
}


# df_results_tmp2_list<-future_apply(df_spectrums_mgf,1, FUN=function(x) search_mgf(x),future.seed=TRUE)
df_results_tmp2_list<-apply(subset,1, search_mgf)

for (i in 1:nrow(subset)){
    z <- search_mgf(subset[i,])
}


df_results_tmp2<-ldply(df_results_tmp2_list,rbind)
new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")

if(output == 1){
    save(df_results_tmp2, file = paste0(new_name,".Rdata"))
}else if(output == 2){
    write.table(df_results_tmp2, file= paste0(new_name,".txt"), row.names=F, col.names=F, quote = FALSE, sep="\t")
}else if(output == 3){
    write_feather(df_results_tmp2, paste0(new_name,".feather"))
}

new<-Sys.time() - old
print(new)
cat("end of search")
