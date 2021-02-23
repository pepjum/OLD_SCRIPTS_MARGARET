
#database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.wavelets_wp.Rdata"
#database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.peaks.Rdata"
database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/human_hcd_selected.peaks_primer_espectro_y_los_que_pasan_filtro_de_rango.Rdata"

#objeto_ruidos<-"/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/objeto_señales_ruido_WP.Rdata"
objeto_ruidos<-"/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/1_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_SPECTRAST.Rdata"

rangeo<-0.1
spread_value<-0
Dalton<-1
output<-1
best_number<-0
#method<-"cos"
#wave<-"wp"

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
#df_spectrums_mgf$idx<-seq(1:nrow(df_spectrums_mgf))
db_spectra<-get(load(database))   #_500-1500.Rdata

plan(multicore, workers=25)

options(future.globals.maxSize= +Inf)

selecting_spectrums_c<-cmpfun(selecting_spectrums)

# names(df_spectrums_mgf_all)[9]<-"Peaks_mgf"
output<-data.frame()
for (i in unique(df_spectrums_mgf_all$idx)){
    df_spectrums_mgf<-df_spectrums_mgf_all[which(df_spectrums_mgf_all$idx== i),]

    print(df_spectrums_mgf$spectrum[1])
    spectra<-df_spectrums_mgf[1,]
    cat("selecting spectrums for MZ range","\n")

    df_spectrums_selected_for_each_spectrum<-ldply(future_apply(spectra,1, function(x) selecting_spectrums_c(x,rangeo,db_spectra), future.seed=TRUE ),rbind)

    # df_spectrums_selected_for_each_spectrum$charge_db<-str_extract(paste(df_spectrums_selected_for_each_spectrum$Name), "[0-9]+")
    # df_spectrums_selected_for_each_spectrum_filtered<-data.frame()
    #
    # for (i in unique(df_spectrums_mgf$idx)){
    #     charge_sel<-df_spectrums_mgf$charge[i]
    #     tmp_selected<-df_spectrums_selected_for_each_spectrum[which(df_spectrums_selected_for_each_spectrum$idx ==i),]
    #     tmp_optimus<-tmp_selected[which(tmp_selected$charge_db == charge_sel),]
    #     df_spectrums_selected_for_each_spectrum_filtered<-rbind(df_spectrums_selected_for_each_spectrum_filtered,tmp_optimus)
    # }

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
    df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(signal,-D)),]
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
    cat("example of results: ",paste(df_results_tmp2[1,]),"\n")
    output<-rbind(output,df_results_tmp2)
    print(dim(output))

}


new_name<-paste(substr(objeto_ruidos, 1, nchar(objeto_ruidos) - 5), "results", sep="")
new_name2<-paste(new_name,method,sep="_")
new_name3<-paste(new_name2,wave,sep="_")

if(output == 1){
    save(output, file = paste0(new_name3,".Rdata"))
}else if(output == 2){
    write.table(output, file= paste0(new_name3,".txt"), row.names=F, col.names=TRUE, quote = FALSE, sep="\t")
}

# output<-df_results_tmp2
# output$signal[1:10]<-paste0(output$signal[1:10],"_",0)
output$noise_level<-unlist(lapply(strsplit(paste(output$signal),"_"),"[",4))
output$noise_level[is.na(output$noise_level)]<-0
output$peaks<-NULL
output$Comments<-NULL
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

processed_output<-data.frame("peptide"=peptide,"position"=position_peptide,"first_hit"=first_hit,"signal"=unique(output$signal),"level"=level)

write.table(processed_output, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/Results_PEAKS_1_peptidos_RANGE_0.1_FULLNIST.txt", col.names=T, row.names=F, sep="\t", quote=F)



peptide_first<-processed_output[which(processed_output$position==1),]

exit<-data.frame()
for (level in unique(peptide_first$level)){
    tmp<-peptide_first[which(peptide_first$level==level),]
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

write.table(exit, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/percentages_of_success_PEAKS_1_peptide_FULL_NIST_0.03.txt", col.names=T, row.names=F, sep="\t", quote=F)

plotter<-unique(exit[,c("level", "percentage")])
names(plotter)[1]<-"sigma"
plotter$SE<-paste("SPECTRAST")
plotter$level<-round(as.numeric(plotter$sigma), digits=4)

library(ggplot2)
pdf("/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/1_peptides_SPECTRAST_good_identified_percentages_by_noise_level_FULL_NIST.pdf")
ggplot(plotter, aes(x = sigma, y =percentage, group=SE, color=SE)) + geom_line() +geom_point()
dev.off()

new<-Sys.time() - old
print(new)
cat("end of search","\n")
