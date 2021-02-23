args=(commandArgs(TRUE))


mgf<-args[1]      #mgf file
database<-args[2]   #database spectaST
rangeo<-args[3]     #for selecting spectra by mz difference of .... (0.1 by default)
spread_value<-args[4]    #spread out
Dalton<-args[5]       # Value of daltons for bining function purposes
output<-args[6]  # 1 for Rdata, 2 for txt, 3 for feather
best_number<-args[7] # select number of matches
method<-args[8]  #dot , dist or cos (dist=euclidean)
wave<- args[9]  # type of wavelet "dwtd", "dwt", "wp"   (donoho+dwt, Discrete wavelet transform, wavelet packets)
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
library(wavethresh)
#setwd(paste0(mgf_dir,"/"))

enableJIT(3)
plan(multicore, workers=10)

options(future.globals.maxSize= +Inf)
cat("reading MGF file","\n")

ReadMGFFile_c<-cmpfun(ReadMGFFile)


mgf_file<-ReadMGFFile_c(mgf)


#READ SPECTRUMS AND KEEP NAMES AND M/Z CHARGE
mgf_to_df_c<-cmpfun(mgf_to_df4)
#df_spectrums_mgf<-mgf_to_df(mgf_file)
df_spectrums_mgf<-mgf_to_df_c(mgf_file)

#df_spectrums_mgf$idx<-seq(1:nrow(df_spectrums_mgf))
rm(mgf_file)

cat("loading database","\n")
db_spectra<-get(load(database))   #_500-1500.Rdata


print("searching for : ", df_spectrums_mgf$spectrum[i])



cat("selecting spectrums for MZ range","\n")

selecting_spectrums_c<-cmpfun(selecting_spectrums_v2)

df_spectrums_selected_for_each_spectrum<-ldply(future_apply(df_spectrums_mgf,1, function(x) selecting_spectrums_c(x,rangeo,db_spectra), future.seed=TRUE ),rbind)

df_spectrums_selected_for_each_spectrum$charge_db<-str_extract(paste(df_spectrums_selected_for_each_spectrum$Name), "[0-9]+")
df_spectrums_selected_for_each_spectrum_filtered<-data.frame()

for (i in unique(df_spectrums_mgf$idx)){
    charge_sel<-df_spectrums_mgf$charge[i]
    tmp_selected<-df_spectrums_selected_for_each_spectrum[which(df_spectrums_selected_for_each_spectrum$idx ==i),]
    tmp_optimus<-tmp_selected[which(tmp_selected$charge_db == charge_sel),]
    df_spectrums_selected_for_each_spectrum_filtered<-rbind(df_spectrums_selected_for_each_spectrum_filtered,tmp_optimus)
}

rm(df_spectrums_selected_for_each_spectrum)
df_spectrums_selected_for_each_spectrum_filtered$.id<-NULL

list_all_peaks<-paste(df_spectrums_mgf$peaks)

#test <- seleccion[, grepl(paste(c("^Name", "PrecursorMZ", "wavelet"), collapse = '|'), colnames(seleccion))]

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


cat("transform to wavelets", "\n")

#wavelets_for_c<-cmpfun(wavelets_for)

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
#base de datos no normalizada
# normalization_w_c<-cmpfun(normalization_w)
# wavelets_norm<-future_lapply(wavelets, FUN = function(x) normalization_w_c(x), future.seed=TRUE)

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$wavelet_mgf<-wavelets

#spectra$idx<-seq(1:(nrow(spectra)))

cat("merge","\n")
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum_filtered, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)
rm(df_spectrums_selected_for_each_spectrum_filtered)

if (method == "dot"){
    cat("calculating dot product and dot bias","\n")

    #save(df_results_tmp2,file="/home/margaret/data/pepe/por_si_fallo_dot_and_db.Rdata")

    Dot_and_DB_comp_w<-cmpfun(Dot_and_DB_waw)

    res<-ldply(future_apply(df_results_tmp2,1, FUN=function(x) Dot_and_DB_comp_w(x), future.seed=TRUE),rbind)

    res<-subset(res, select= -.id)

    df_results_tmp2<-cbind(df_results_tmp2,res)

    df_results_tmp2$wavelet_mgf<-NULL
    df_results_tmp2$wavelet<- NULL
    df_results_tmp2<-df_results_tmp2[which(as.numeric(df_results_tmp2$D) > 0),]
    df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,-D)),]

    ###select only the 10 best matches
    #best_number<-10

    if (best_number !=0){
        cat("selecting only the best  ", best_number, "PSM","\n")
        selected_data<-list()
        for (i in unique(df_results_tmp2$idx)){
            tmp<-df_results_tmp2[which(df_results_tmp2$idx == i),]
            if(nrow(tmp)< best_number){
                tmp_selected<-tmp
                tmp_selected <- na.omit(tmp_selected)
            }else{
                tmp_selected<-tmp[1:best_number,]
                tmp_selected <- na.omit(tmp_selected)
            }
            selected_data[[i]]<-tmp_selected
        }
        df_results_tmp2<-ldply(selected_data,rbind)
    }

    cat("calculating DeltaD","\n")
    #df_results_tmp2$Delta<-NA

    df_results_tmp2[which(is.nan(df_results_tmp2$DB)), "DB"]<-0

    Delta_c<-cmpfun(Delta)

    df_results_tmp3<-Delta_c(df_results_tmp2)
    rm(df_results_tmp2)

    cat("calculating Fscore","\n")

    F_score_com<-cmpfun(F_score)

    res_2<-future_apply(df_results_tmp3,1, FUN = function(x) F_score(x), future.seed=TRUE)
    res_2<-as.data.frame(res_2)
    names(res_2)<-"Fscore"

    df_results_tmp2<-cbind(df_results_tmp3,res_2)

    new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), df_spectrums_mgf$spectrum[i], sep="")
    new_name2<-paste(new_name, "dot", sep="_")
    new_name3<-paste(new_name2,wave,sep="_")

    if(output == 1){
        save(df_results_tmp2, file = paste0(new_name3,".Rdata"))
    }else if(output == 2){
        write.table(df_results_tmp2, file= paste0(new_name3,".txt"), row.names=F, col.names=F, quote = FALSE, sep="\t")
    }

}else if(method=="dist"){
    dist_wav_c<-cmpfun(Dist_wav)
    dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x), future.seed=TRUE)
    df_results_tmp2$Dist<-dist_list
    df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,Dist)),]
    df_results_tmp2$wavelet<-NULL
    df_results_tmp2$wavelet_mgf<-NULL


    if (best_number !=0){
        cat("selecting only the best ", best_number, "PSM for each spectrum","\n")
        selected_data<-list()
        for (i in unique(df_results_tmp2$idx)){
            tmp<-df_results_tmp2[which(df_results_tmp2$idx == i),]
            tmp_selected<-tmp[1:best_number,]
            tmp_selected <- na.omit(tmp_selected)
            selected_data[[i]]<-tmp_selected
        }
        df_results_tmp2<-ldply(selected_data,rbind)
    }


    new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), df_spectrums_mgf$spectrum[i], sep="_")
    new_name2<-paste(new_name,"dist",sep="_")
    new_name3<-paste(new_name2,wave,sep="_")

    if(output == 1){
        save(df_results_tmp2, file = paste0(new_name3,".Rdata"))
    }else if(output == 2){
        write.table(df_results_tmp2, file= paste0(new_name3,".txt"), row.names=F, col.names=TRUE, quote = FALSE, sep="\t")
    }

}else if(method=="cos"){
    library(stylo)
    Dist_cosine_wav_c<-cmpfun(Dist_cosine_wav_v2)
    dist_cosine<-future_apply(df_results_tmp2,1, FUN = function(x) Dist_cosine_wav_c(x, db_spectra), future.seed=TRUE)
    df_results_tmp2$Dist_cosines<-paste(dist_cosine)
    df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,Dist_cosines)),]
    df_results_tmp2$wavelet<-NULL
    df_results_tmp2$wavelet_mgf<-NULL


    if (best_number !=0){
        cat("selecting only the best ", best_number, "PSM for each spectrum","\n")
        selected_data<-list()
        for (i in unique(df_results_tmp2$idx)){
            tmp<-df_results_tmp2[which(df_results_tmp2$idx == i),]
            tmp_selected<-tmp[1:best_number,]
            tmp_selected <- na.omit(tmp_selected)
            selected_data[[i]]<-tmp_selected
        }
        df_results_tmp2<-ldply(selected_data,rbind)
    }


    new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), df_spectrums_mgf$idx[i], sep="")
    new_name2<-paste(new_name,"cos",sep="_")
    new_name3<-paste(new_name2,wave,sep="_")

    if(output == 1){
        save(df_results_tmp2, file = paste0(new_name3,".Rdata"))
    }else if(output == 2){
        write.table(df_results_tmp2, file= paste0(new_name3,".txt"), row.names=F, col.names=F, quote = FALSE, sep="\t")
    }

}else if(method=="mutual"){
    library(entropy)
    mut_wav_c<-cmpfun(mutual_information)
    mutual_inf<-future_apply(df_results_tmp2,1, FUN = function(x) mut_wav_c(x), future.seed=TRUE)

    df_results_tmp2$mutual_inf<-paste(mutual_inf)
    df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,mutual_inf)),]
    df_results_tmp2$wavelet<-NULL
    df_results_tmp2$wavelet_mgf<-NULL

    if (best_number !=0){
        cat("selecting only the best ", best_number, "PSM for each spectrum","\n")
        selected_data<-list()
        for (i in unique(df_results_tmp2$idx)){
            tmp<-df_results_tmp2[which(df_results_tmp2$idx == i),]
            tmp_selected<-tmp[1:best_number,]
            tmp_selected <- na.omit(tmp_selected)
            selected_data[[i]]<-tmp_selected
        }
        df_results_tmp2<-ldply(selected_data,rbind)
    }


    new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")
    new_name2<-paste(new_name,"mutual_inf",sep="_")
    new_name3<-paste(new_name2,wave,sep="_")

    if(output == 1){
        save(df_results_tmp2, file = paste0(new_name3,".Rdata"))
    }else if(output == 2){
        write.table(df_results_tmp2, file= paste0(new_name3,".txt"), row.names=F, col.names=F, quote = FALSE, sep="\t")
    }


}


new<-Sys.time() - old
print(new)
cat("end of search","\n")
