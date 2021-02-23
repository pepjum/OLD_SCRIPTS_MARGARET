source("/home/margaret/data/pepe/scripts/functions_SE.R")
#mgf_dir="."
cat("started at","\n")
old<-Sys.time()
print(old)

options(stringsAsFactors = FALSE)

rangeo<-0.1 # cambiar a range =1
spread_value<-0
Dalton<-1
n_peaks<-6
intensity<-500

library(stringr)
library(stringi)
library(plyr)
library(future.apply)
library(compiler)
library(wavethresh)
#setwd(paste0(mgf_dir,"/"))

enableJIT(3)
plan(multicore, workers=40)

options(future.globals.maxSize= +Inf)
cat("reading MGF file","\n")

objeto_seleccion<-get(load("~/data/pepe/SEARCH_ENGINE_DIC20/100_PEPTIDES/100_peptidos_seleccionados_DB.Rdata"))

cat("creating a vector adding zeros to positions where we have not mz values","\n")

picos<-paste(objeto_seleccion$Peaks)

descomprimir<-function(df_results_tmp2){
    tmp_mgf<-strsplit(paste(df_results_tmp2), "###")
    m_mgf<-paste(lapply(strsplit(unlist(tmp_mgf),"\t"),'[',1))
    z_mgf<-paste(lapply(strsplit(unlist(tmp_mgf),"\t"),'[',2))
    df_peaks_mgf<-data.frame("M"=as.numeric(m_mgf),"normalized"=as.numeric(z_mgf))
}


list_df_matrix<-future_lapply(picos, FUN = function(x) descomprimir(x), future.seed=TRUE)

original<-(list_df_matrix)

end_SNR_Tot<-list()
end_PSI_signals_Tot<-list()
end_PSI_noises_Tot<-list()
end_signals_Tot<-list()

for(l in 1:length(original)){
    cat("peptido",l,"\n")
    mu<-0
    sigma<-seq(0.01,0.3, length.out=15)

    SNR_signals<-list()
    Psi_signals<-list()
    Psi_noises<<-list()
    signals<-list()
    signal<-original[l]
    signal<-as.data.frame(signal)
    for (i in sigma){
        print(i)
        for(j in seq(1,15)){
            noise<-rnorm(nrow(signal), mu,i)
            noise[which(noise <0)]<-0
            noise[which(noise >1)]<-1
            Psi_noise<-sum(noise^2)
            Psi_noises[[paste0('peptide_',l,'_Noise_',i, '_',j)]]<-Psi_noise
            signal$normalized<-as.numeric(signal$normalized) + noise
            signal$normalized[which(signal$normalized >1)]<-1
            Psi_signal<-sum(signal$normalized^2)
            Psi_signals[[paste0('peptide_', l,'_Noise_',i, '_',j)]]<-Psi_signal
            SNR<-Psi_signal/Psi_noise
            SNR_signals[[paste0('peptide_', l,'_Noise_',i, '_',j)]]<-SNR
            signals[[paste0('peptide_', l,'_Noise_',i, '_',j)]] <- signal
        }
        end_SNR_tmp<-list()
        end_SNR_tmp[[1]]<-Inf
        end_SNR_signal<-c(end_SNR_tmp,SNR_signals)

        end_PSI_signal_tmp<-list()
        end_PSI_signal_tmp[[1]]<-1
        end_PSI_signal<-c(end_PSI_signal_tmp,Psi_signals)

        end_PSI_noise_tmp<-list()
        end_PSI_noise_tmp[[1]]<-0
        end_PSI_noise<-c(end_PSI_noise_tmp,Psi_noises)

        end_list<-list()
        end_list[[1]]<-original[l]
        end_list<-c(end_list,signals)

    }
    end_SNR_Tot<-c(end_SNR_Tot, end_SNR_signal)
    end_PSI_signals_Tot<-c(end_PSI_signals_Tot, end_PSI_signal)
    end_PSI_noises_Tot<-c(end_PSI_noises_Tot, end_PSI_noise)
    end_signals_Tot<-c(end_signals_Tot, end_list)
}


save(end_signals_Tot,file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/100_peptidos_señales_con_ruido_mu_0_sigma_0.1_to_0.3.Rdata")

# adapto el dataframe para poder tener la misma estructura que en los SE programados
df_spectrums_mgf_bck<-objeto_seleccion
# me quedo solo con el peptido al que le vamos a añadir ruido
#df_spectrums_mgf<-df_spectrums_mgf[2,]
df_spectrums_mgf$Peaks<-NULL
df_spectrums_mgf$idx<-seq(1:100)
names(df_spectrums_mgf)[1:2]<-c("spectrum","mz_value")

# repito los datos de identificacion del peptido el mismo numero de veces que señales hemos simulado
df_spectrums_mgf<-df_spectrums_mgf[rep(seq_len(nrow(df_spectrums_mgf)), each=226),]   #226=15*15+1
df_spectrums_mgf$signal<-names(end_signals_Tot)
df_spectrums_mgf[which(df_spectrums_mgf$signal==""),"signal"]<-paste0("peptide_",seq(1:100))
#guardo la lista original por si acaso
list_df_matrix_bck<-list_df_matrix
#nombro como el objeto de siempre la lista de las señales con ruido para seguir trabajando con el codigo viejo

list_df_matrix<-end_signals_Tot

# señal para SPECTRAST

cat("convert peaks dataframe into string to indexing in the info dataframe","\n")

codification_c<-cmpfun(codification)

list_linecoded<-future_lapply(list_df_matrix, FUN = function(x) codif100_peptidos_señales_con_ruido_mu_0_sigma_0ication_c(x), future.seed=TRUE)

df_spectrums_mgf$PSI_signal<-paste(end_PSI_signals_Tot)
df_spectrums_mgf$PSI_noise<-paste(end_PSI_noises_Tot)
df_spectrums_mgf$SNR<-paste(end_SNR_Tot)
df_spectrums_mgf$Peaks_mgf<-paste(list_linecoded)

save(df_spectrums_mgf, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/100_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_SPECTRAST.Rdata")

df_spectrums_mgf$Peaks_mgf<-NULL

# señal para WP

cat("transform to wavelets", "\n")

divide<-list()
counter<-1
for (index in seq(1,length(list_df_matrix),(length(list_df_matrix)/100))){
    print(index)
    if(index+30009 <length(list_df_matrix)){
        lista<-list_df_matrix[index:(index+30009)]
        divide[[counter]]<-lista
        counter<-counter+1

    }else{
        lista<-list_df_matrix[index:length(list_df_matrix)]
        divide[[counter]]<-lista
        counter<-counter+1
    }
}

wavelets_c<-cmpfun(wavelets)

# output_wavelets<-list()

# for(i in 1:length(divide)){
    # cat("section", i, "\n")
    # list_df_matrix<-divide[[i]]
    wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)
    # output_wavelets<-c(output_wavelets,wavelets)
# }

wavelets<-output_wavelets



df_spectrums_mgf$wavelet_mgf<-wavelets

save(df_spectrums_mgf, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/100_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_WP.Rdata")


wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)




#señal para WD
df_spectrums_mgf$wavelet_mgf<-NULL

cat("transform to wavelets", "\n")

# divide<-list()
# counter<-1
# for (index in seq(1,length(list_df_matrix)-10000,10001)){
#     print(index)
#     lista<-list_df_matrix[index:(index+10000)]
#     divide[[counter]]<-lista
#     counter<-counter+1
# }

wavelets_c<-cmpfun(wavelets_d)
# output_wavelets<-list()

# for(i in 1:length(divide)){
    # cat("section", i, "\n")
    # list_df_matrix<-divide[[i]]
    wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)
    # output_wavelets<-c(output_wavelets,wavelets)
# }

# wavelets<-output_wavelets

df_spectrums_mgf$wavelet_mgf<-wavelets

save(df_spectrums_mgf, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/100_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_DWT.Rdata")




wavelets_c<-cmpfun(wavelets_d)
# output_wavelets<-list()

# for(i in 1:length(divide)){
    # cat("section", i, "\n")
    # list_df_matrix<-divide[[i]]
    wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)
    # output_wavelets<-c(output_wavelets,wavelets)
# }

# wavelets<-output_wavelets

df_spectrums_mgf$wavelet_mgf<-wavelets

save(df_spectrums_mgf, file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/100_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_DWT.Rdata")



wavelets_c<-cmpfun(wavelets_for)

wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)


df_spectrums_mgf$wavelet_mgf<-wavelets

save(df_spectrums_mgf, file="/mnt/md0/pepe/SEARCH_ENGINE_DIC20/100_PEPTIDES/100_peptidos_objeto_señales_ruido__sigma0.01_to_0.3_DWTD.Rdata")
