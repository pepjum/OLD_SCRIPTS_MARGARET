#database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.wavelets_dwt.Rdata"
database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.wavelets_wp.Rdata"

#objeto_ruidos<-"/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/objeto_señales_ruido_WP.Rdata"
#objeto_ruidos<-"/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/100_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_DWT.Rdata"
objeto_ruidos<-"/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/100_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_WP.Rdata"


rangeo<-0.05
spread_value<-0
Dalton<-1
output<-1
best_number<-0
method<-"cos"
#wave<-"dwt"
wave<-"wp"

source("/home/margaret/data/pepe/scripts/functions_SE.R")

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
plan(multicore, workers=30)

options(future.globals.maxSize= +Inf)

df_spectrums_mgf_all<-get(load(objeto_ruidos))
db_spectra<-get(load(database))   #_500-1500.Rdata
#df_spectrums_mgf_all$peaks<-NULL

output<-data.frame()

for (i in unique(df_spectrums_mgf_all$idx)){
    df_spectrums_mgf<-df_spectrums_mgf_all[which(df_spectrums_mgf_all$idx== i),]
    print(df_spectrums_mgf$spectrum[1])
    spectra<-df_spectrums_mgf[1,]
    cat("selecting spectrums for MZ range","\n")

    selecting_spectrums_c<-cmpfun(selecting_spectrums_v2)

    df_spectrums_selected_for_each_spectrum<-ldply(future_apply(spectra,1, function(x) selecting_spectrums_c(x,rangeo,db_spectra), future.seed=TRUE ),rbind)
    cat("peptides in DB for peptide ",i,": ", nrow(df_spectrums_selected_for_each_spectrum),"\n")
    #df_spectrums_selected_for_each_spectrum$charge_db<-str_extract(paste(df_spectrums_selected_for_each_spectrum$Name), "[0-9]+")

    cat("merge","\n")
    df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )
    cat("merged: ",nrow(df_results_tmp),"\n" )
    #eliminar no hits
    df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
    rm(df_results_tmp)
    rm(df_spectrums_selected_for_each_spectrum_filtered)

    if (method == "dot"){
        cat("calculating dot product and dot bias","\n")

        #save(df_results_tmp2,file="/home/margaret/data/pepe/por_si_fallo_dot_and_db.Rdata")

        Dot_and_DB_comp_w<-cmpfun(Dot_and_DB)

        res<-ldply(future_apply(df_results_tmp2,1, FUN=function(x) Dot_and_DB_comp_w(x, db_spectra), future.seed=TRUE),rbind)

        res<-subset(res, select= -.id)

        df_results_tmp2<-cbind(df_results_tmp2,res)

        df_results_tmp2$wavelet_mgf<-NULL

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

    output<-rbind(output,df_results_tmp2)
    print(dim(output))
    }else if(method=="dist"){
        dist_wav_c<-cmpfun(Dist_wav_v2)
        dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x, db_spectra), future.seed=TRUE)
        df_results_tmp2$Dist<-dist_list
        df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,as.numeric(Dist))),]
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

    output<-rbind(output, df_results_tmp2)
    print(dim(output))

    }else if(method=="cos"){
        library(stylo)
        cat("calculating distances", "\n")
        Dist_cosine_wav_c<-cmpfun(Dist_cosine_wav_v2)
        dist_cosine<-future_apply(df_results_tmp2,1, FUN = function(x) Dist_cosine_wav_c(x, db_spectra), future.seed=TRUE)
        df_results_tmp2$Dist_cosines<-paste(dist_cosine)
        df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(signal,as.numeric(Dist_cosines))),]
        # df_results_tmp2$wavelet<-NULL
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

    output<-rbind(output, df_results_tmp2)
    print(dim(output))
    cat("example of results:", paste(df_results_tmp2[1,]),"\n")

    }else if(method=="mutual"){
        library(entropy)
        mut_wav_c<-cmpfun(mutual_information_v2)
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

    output<-rbind(output, df_results_tmp2)

    }
}
