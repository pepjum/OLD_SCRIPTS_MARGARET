
#database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.wavelets_wp.Rdata"
database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/human_hcd_selected._wavelets_PACKETS_one_spectra_and_the_filtered_by_range_wp_full.Rdata"

#objeto_ruidos<-"/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/objeto_señales_ruido_WP.Rdata"
objeto_ruidos<-"/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/1_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_WP.Rdata"


rangeo<-0.1
spread_value<-0
Dalton<-1
output<-1
best_number<-0
method<-"cos"
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

output<-data.frame()

for (i in unique(df_spectrums_mgf_all$idx)){

    df_spectrums_mgf<-df_spectrums_mgf_all[which(df_spectrums_mgf_all$idx== i),]


    #print(i)
    spectra<-df_spectrums_mgf[i,]
    #spectra$idx<-NULL
    print(df_spectrums_mgf$spectrum[1])

    cat("selecting spectrums for MZ range","\n")

    selecting_spectrums_c<-cmpfun(selecting_spectrums_v2)

    df_spectrums_selected_for_each_spectrum<-ldply(future_apply(spectra,1, function(x) selecting_spectrums_c(x,rangeo,db_spectra), future.seed=TRUE ),rbind)

    df_spectrums_selected_for_each_spectrum$charge_db<-str_extract(paste(df_spectrums_selected_for_each_spectrum$Name), "[0-9]+")
    # df_spectrums_selected_for_each_spectrum_filtered<-data.frame()
    cat("spectrums in dB: ", nrow(df_spectrums_selected_for_each_spectrum),"\n")
    # for (i in unique(df_spectrums_mgf$idx)){
    #     charge_sel<-df_spectrums_mgf$charge[i]
    #     tmp_selected<-df_spectrums_selected_for_each_spectrum[which(df_spectrums_selected_for_each_spectrum$idx ==i),]
    #     tmp_optimus<-tmp_selected[which(tmp_selected$charge_db == charge_sel),]
    #     df_spectrums_selected_for_each_spectrum_filtered<-rbind(df_spectrums_selected_for_each_spectrum_filtered,tmp_optimus)
    # }

    # rm(df_spectrums_selected_for_each_spectrum)


    cat("merge","\n")
    df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )

    #eliminar no hits
    df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
    rm(df_results_tmp)
    #rm(df_spectrums_selected_for_each_spectrum_filtered)

    if (method == "dot"){
        cat("calculating dot product and dot bias","\n")

        #save(df_results_tmp2,file="/home/margaret/data/pepe/por_si_fallo_dot_and_db.Rdata")

        Dot_and_DB_comp_w<-cmpfun(Dot_and_DB)

        res<-ldply(future_apply(df_results_tmp2,1, FUN=function(x) Dot_and_DB_comp_w(x, db_spectra), future.seed=TRUE),rbind)

        res<-subset(res, select= -.id)

        df_results_tmp2<-cbind(df_results_tmp2,res)

        df_results_tmp2$wavelet_mgf<-NULL

        df_results_tmp2<-df_results_tmp2[which(as.numeric(df_results_tmp2$D) > 0),]
        df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(signal,-D)),]

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

    }else if(method=="dist"){
        dist_wav_c<-cmpfun(Dist_wav_v2)
        dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x, db_spectra), future.seed=TRUE)
        df_results_tmp2$Dist<-dist_list
        df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(signal,as.numeric(Dist))),]
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


new_name<-paste(substr(objeto_ruidos, 1, nchar(objeto_ruidos) - 5), "results", sep="")
new_name2<-paste(new_name,method,sep="_")
new_name3<-paste(new_name2,wave,sep="_")

if(output == 1){
    save(output, file = paste0(new_name3,".Rdata"))
}else if(output == 2){
    write.table(output, file= paste0(new_name3,".txt"), row.names=F, col.names=TRUE, quote = FALSE, sep="\t")
}

output$noise_level<-unlist(lapply(strsplit(paste(output$signal),"_"),"[",4))
output$noise_level[is.na(output$noise_level)]<-0


peptide<-c()
position_peptide<-c()
first_hit<-c()
level<-c()

for(signal in unique(output$signal)){
    print(signal)
    tmp<-output[which(output$signal ==signal),]
    first_hit<-c(first_hit,paste(tmp$Name)[1])
    for_position<-paste(tmp$spectrum[1])
    peptide<-c(peptide, for_position)
    position<-which(tmp$Name==for_position)
    if (length(position)==0){
        position<-"NA"
    }
    position_peptide<-c(position_peptide,position)
    level<-c(level, tmp$noise_level[[1]][1])
}

processed_output<-data.frame("peptide"=peptide,"position"=position_peptide,"first_hit"=first_hit,"signal"=unique(output$signal),"sigma"=level)

write.table(processed_output, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/WP_1_peptidos_RANGE_0.1_FULL_NIST.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(processed_output, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/DWT_100_peptidos_RANGE_0.05.txt", col.names=T, row.names=F, sep="\t", quote=F)


peptide_first<-processed_output[which(processed_output$position==1),]

exit<-data.frame()
for (level in unique(peptide_first$sigma)){
    tmp<-peptide_first[which(peptide_first$sigma==level),]
    count<-nrow(tmp)
    tmp$count<-count
    if(level==0){
        percentage<-(count*100)/1
        tmp$percentage<-percentage
    }else{
        percentage<-(count*100)/15
        tmp$percentage<-percentage
    }
    exit<-rbind(exit,tmp)
}

write.table(exit, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/percentages_of_success_WP_100_peptides_0.03.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(exit, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/percentages_of_success_DWT_100_peptides_0.03.txt", col.names=T, row.names=F, sep="\t", quote=F)

plotter_DWT<-unique(exit_DWT[,c("level", "percentage")])
plotter_DWT$SE<-paste("DWT")
plotter_DWT$level<-round(as.numeric(plotter_DWT$level), digits=4)


plotter_WP<-unique(exit[,c("level", "percentage")])
plotter_WP$SE<-paste("WP")
plotter_WP$level<-round(as.numeric(plotter_WP$level), digits=4)

plotter_tot<-rbind(plotter_WP, plotter_DWT)
library(ggplot2)
pdf("/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/100_peptides_DWT_&_WP_good_identified_percentages_by_noise_level.pdf")
ggplot(plotter_tot, aes(x = level, y =percentage, group=SE, color=SE)) + geom_line() +geom_point()
dev.off()



new<-Sys.time() - old
print(new)
cat("end of search","\n")
