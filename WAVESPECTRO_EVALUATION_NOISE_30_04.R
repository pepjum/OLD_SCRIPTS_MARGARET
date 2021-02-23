source("/home/margaret/data/pepe/scripts/functions_SE.R")
#mgf_dir="."
cat("started at","\n")
old<-Sys.time()
print(old)

options(stringsAsFactors = FALSE)

mgf<-"/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/MS2PIP_Predictions.mgf"

#database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.wavelets_packets.Rdata"
database<-"/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.wavelets_wp.Rdata"


rangeo<-0.01 # cambiar a range =1
spread_value<-0
Dalton<-1
n_peaks<-150
intensity<-500

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

rm(mgf_file)

cat("loading database","\n")
db_spectra<-get(load(database))   #_150-1500.Rdata

cat("selecting spectrums for MZ range","\n")

selecting_spectrums_c<-cmpfun(selecting_spectrums)

df_spectrums_selected_for_each_spectrum<-ldply(future_apply(df_spectrums_mgf,1, function(x) selecting_spectrums_c(x,rangeo,db_spectra), future.seed=TRUE ),rbind)

library(stringr)

df_spectrums_selected_for_each_spectrum$charge_db<-str_extract(paste(df_spectrums_selected_for_each_spectrum$Name), "[0-9]+")
df_spectrums_selected_for_each_spectrum_filtered<-data.frame()

for (i in unique(df_spectrums_mgf$idx)){
    charge_sel<-df_spectrums_mgf$charge[i]
    tmp_selected<-df_spectrums_selected_for_each_spectrum[which(df_spectrums_selected_for_each_spectrum$idx ==i),]
    tmp_optimus<-tmp_selected[which(tmp_selected$charge_db == charge_sel),]
    df_spectrums_selected_for_each_spectrum_filtered<-rbind(df_spectrums_selected_for_each_spectrum_filtered,tmp_optimus)
}

list_all_peaks<-paste(df_spectrums_mgf$peaks)

#test <- seleccion[, grepl(paste(c("^Name", "PrecursorMZ", "wavelet"), collapse = '|'), colnames(seleccion))]

cat("transform peaks info to a dataframe to work in it","\n")

convert_df_c<-cmpfun(convert_df)

list_df_peaks_not_binned<-future_lapply(list_all_peaks, FUN = function(x) convert_df_c(x), future.seed=TRUE)

### sacar señal en plot
library(ggplot2)

señal_raw<-as.data.frame(list_df_peaks_not_binned[2])
señal_raw$signal<-rep("raw",nrow(señal_raw))
señal_raw$M<-as.numeric(señal_raw$M)
señal_raw$Z<-as.numeric(señal_raw$Z)
pdf("signal_raw.pdf")
ggplot(data=señal_raw, aes(x=M, y=Z, fill=signal)) + geom_line()
dev.off()


cat("bining","\n")

bining_c<-cmpfun(bining)
list_df_peaks_binned<-future_lapply(list_df_peaks_not_binned, FUN= function(x) bining_c(x,1), future.seed=TRUE)

señal_binada<-as.data.frame(list_df_peaks_binned[2])
señal_binada$signal<-rep("binada",nrow(señal_binada))
señal_binada$M<-as.numeric(señal_binada$M)
señal_binada$Z<-as.numeric(señal_binada$Z)
names(señal_raw)<-c("M","Intensity","signal")
names(señal_binada)<-c("M","Intensity","signal")

plotter<-data.frame("M"<-NULL,"Intensity"<-NULL,"signal"<-NULL)
plotter<-rbind(plotter,señal_raw)
plotter<-rbind(plotter,señal_binada)

pdf("señales_simuladas.pdf")
ggplot(data=plotter, aes(x=M, y=Intensity, fill=signal, color=signal)) + geom_line()
dev.off()


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

señal_mgf_rango<-as.data.frame(list_selected[2])
names(señal_mgf_rango)[2]<-"Intensity"
señal_mgf_rango$signal<-rep("LGFLHSGTAK",nrow(señal_mgf_rango))



normalization_c<-cmpfun(normalization)

list_df_normalized<-future_lapply(list_selected, FUN = function(x) normalization_c(x), future.seed=TRUE)

señal_mgf_normalizada<-as.data.frame(list_df_normalized[2])
names(señal_mgf_normalizada[2])<-"Intensity"
señal_mgf_normalizada$signal<-rep("LGFLHSGTAK_norm",nrow(señal_mgf_normalizada))

plotter<-data.frame("M"<-NULL,"Intensity"<-NULL,"signal"<-NULL)
plotter<-rbind(plotter,señal_mgf_rango)
plotter<-rbind(plotter,señal_mgf_normalizada)

pdf("señal_LGFLHSGTAK_selected_and_normalized.pdf")
ggplot(data=plotter, aes(x=M, y=Intensity, fill=signal, color=signal)) + geom_line()
dev.off()


cat("creating a vector adding zeros to positions where we have not mz values","\n")

create_full_vector_c<-cmpfun(create_full_vector)

list_df_matrix<-future_lapply(list_df_normalized, FUN = function(x) create_full_vector_c(x), future.seed=TRUE)

vector_signal<-as.data.frame(list_df_matrix[2])
names(vector_signal)[2]<-"Intensity"
vector_signal$señal<-rep("vectorizada",nrow(vector_signal))

pdf("signal_vectorized.pdf")
ggplot(data=vector_signal, aes(x=M, y=Intensity, fill=señal)) + geom_line()
dev.off()


########## hasta aqui lo comun
### cargar la base de datos de picos para extraer el espectrio LGFLHSGTAK y el de mayor hit en resultados de picos

database<-get(load("/home/margaret/data/pepe/HUMAN_NIST_25_7_2018/HUMAN_best.peaks.Rdata"))

espectro_LGFLHSGTAK<-database[which(database$Name=="LGFLHSGTAK/2"),]

resultados<-get(load("/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/MS2PIP_Predictions.results_dist_wp_original_wp.Rdata"))
######### parte de wavelets
espectro_LTAEELLNK<-database[which(database$Name=="LTAEELLNK/2"),]
espectro_LTAEELLNK<-database[which(database$Name=="LTAEELLNK/2"),]
tmp_database_e1<-strsplit(paste(espectro_LGFLHSGTAK["Peaks"]), "###")
m_database<-paste(lapply(strsplit(unlist(tmp_database_e1),"\t"),'[',1))
z_database<-paste(lapply(strsplit(unlist(tmp_database_e1),"\t"),'[',2))
df_peaks_database_LGFLHSGTAK<-data.frame("M"=as.numeric(m_database),"Z"=as.numeric(z_database))

tmp_database_e2<-strsplit(paste(espectro_LTAEELLNK["Peaks"]), "###")
m_database<-paste(lapply(strsplit(unlist(tmp_database_e2),"\t"),'[',1))
z_database<-paste(lapply(strsplit(unlist(tmp_database_e2),"\t"),'[',2))
df_peaks_database_LTAEELLNK<-data.frame("M"=as.numeric(m_database),"Z"=as.numeric(z_database))

tmp_database_e3<-strsplit(paste(espectro_LTAEELLNK["Peaks"]), "###")
m_database<-paste(lapply(strsplit(unlist(tmp_database_e3),"\t"),'[',1))
z_database<-paste(lapply(strsplit(unlist(tmp_database_e3),"\t"),'[',2))
df_peaks_database_LTAEELLNK<-data.frame("M"=as.numeric(m_database),"Z"=as.numeric(z_database))

df_peaks_database_LGFLHSGTAK$signal<-rep("LGFLHSGTAK_DB",nrow(df_peaks_database_LGFLHSGTAK))
names(df_peaks_database_LGFLHSGTAK)[2]<-"Intensity"
df_peaks_database_LTAEELLNK$signal<-rep("LTAEELLNK_DB",nrow(df_peaks_database_LTAEELLNK))
names(df_peaks_database_LTAEELLNK)[2]<-"Intensity"
df_peaks_database_LTAEELLNK$signal<-rep("LTAEELLNK_DB",nrow(df_peaks_database_LTAEELLNK))
names(df_peaks_database_LTAEELLNK)[2]<-"Intensity"



pdf("signal_simulated+database+most_similar.pdf")
par(mfrow=c(3,1))

ggplot(data=vector_signal, aes(x=M, y=Intensity, fill=señal)) + geom_line()
ggplot(data=df_peaks_database_LGFLHSGTAK, aes(x=M, y=Intensity, fill=signal)) + geom_line()
ggplot(data=df_peaks_database_LTAEELLNK, aes(x=M, y=Intensity, fill=signal)) + geom_line()

dev.off()

vector_signal$Intensity<-(vector_signal$Intensity)*(-1)
names(vector_signal)[3]<-"signal"
#df_peaks_database_LTAEELLNK$Intensity<-(df_peaks_database_LTAEELLNK$Intensity)*(-1)
df_peaks_database_LTAEELLNK_inv<-df_peaks_database_LTAEELLNK
df_peaks_database_LTAEELLNK_inv$Intensity<-(df_peaks_database_LTAEELLNK_inv$Intensity)*(-1)
plotter_1<-data.frame()
plotter_2<-data.frame()
plotter_3<-data.frame()
plotter_1<-rbind(vector_signal,df_peaks_database_LGFLHSGTAK)
plotter_2<-rbind(vector_signal,df_peaks_database_LTAEELLNK)
plotter_3<-rbind(df_peaks_database_LGFLHSGTAK,df_peaks_database_LTAEELLNK_inv)

pdf("signals_plotted_comparisons_dist_cosines_wp_good_identification.pdf")
#par(mfrow=c(3,1))
ggplot(data=plotter_1, aes(x=M, y=Intensity, fill=signal, color=signal)) + geom_line()
ggplot(data=plotter_2, aes(x=M, y=Intensity, fill=signal, color=signal)) + geom_line()
ggplot(data=plotter_3, aes(x=M, y=Intensity, fill=signal, color=signal)) + geom_line()

dev.off()
cat("transform to wavelets", "\n")

#wavelets_c<-cmpfun(wavelets) #esto es para wp
wavelets_c<-cmpfun(wavelets) #esto es para wd

wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$wavelet_mgf<-wavelets

df_spectrums_mgf$idx<-seq(1:(nrow(df_spectrums_mgf)))

cat("merge","\n")
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum_filtered, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)

### saving object
save(df_results_tmp2,file="/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/selection_TP53_ms2pip_after_distances_wp.Rdata")
###### cambiar para wavelets_packets
cat("calculating distances betweeen cosines", "\n")

dist_wav_c<-cmpfun(Dist_cosine_wav)
dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x), future.seed=TRUE)
df_results_tmp2$Dist<-as.numeric(paste(dist_list))
df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,Dist)),]
df_results_tmp2$wavelet<-NULL
df_results_tmp2$wavelet_mgf<-NULL

best_number<-0
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

cat("saving results", "\n")
wave<-"wp"
new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")
new_name2<-paste(new_name,"dist",sep="_")
new_name3<-paste(new_name2,wave,sep="_")
new_name4<-paste0(new_name3,"_original_wp")

save(df_results_tmp2, file = paste0(new_name4,".Rdata"))

##### informacion mutual
library(entropy)

mut_wav_c<-cmpfun(mutual_information)

mutual_inf<-future_apply(df_results_tmp2,1, FUN = function(x) mut_wav_c(x), future.seed=TRUE)

df_results_tmp2$mutual_inf<-as.numeric(paste(mutual_inf))
df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,mutual_inf)),]
df_results_tmp2$wavelet<-NULL
df_results_tmp2$wavelet_mgf<-NULL




#############################################################################################################
###############################################################################################################
######## calculate signal to noise ratio in original simulated signal
#############################################################################################################3
##############################################################################################################

#this is the signal
library("MassSpecWavelet")

cat("reading MGF file","\n")

ReadMGFFile_c<-cmpfun(ReadMGFFileMASSSimulator)


mgf_file<-ReadMGFFile_c(mgf)


#READ SPECTRUMS AND KEEP NAMES AND M/Z CHARGE
mgf_to_df_c<-cmpfun(mgf_to_df2)
#df_spectrums_mgf<-mgf_to_df(mgf_file)
df_spectrums_mgf<-mgf_to_df_c(mgf_file)

list_all_peaks<-paste(df_spectrums_mgf$peaks)

tmp<-strsplit(paste(list_all_peaks), "###")
m<-paste(lapply(strsplit(unlist(tmp)," "),'[',1))
z<-paste(lapply(strsplit(unlist(tmp)," "),'[',2))
df<-data.frame("z"=z)
matriz<-data.matrix(df)

#DOING CWT
results_to_plot<-data.frame("sum(peakSNR)"=NULL,"dist_cosines"=NULL, "signal"<-NULL)

#data(exampleMS)
SNR.Th <- 2
peakInfo <- peakDetectionCWT(matriz,scales=1:30, SNR.Th=SNR.Th)
majorPeakInfo = peakInfo$majorPeakInfo
peakIndex <- majorPeakInfo$peakIndex

list0<-paste(peakInfo$majorPeakInfo$peakSNR)   #este identifica
distancia_id0<-as.numeric(df_results_tmp2$Dist[1])  #3.29

original<-data.frame("sum(peakSNR)"=sum(as.numeric(list0)),"dist_cosines"=distancia_id0, "signal"= 0)

results_to_plot<-rbind(results_to_plot, original)
##########################################################################################################
##########################################################################################################l
###### ESPECTRO añadiendo más ruido
library(RMThreshold)

#la señal original está en el objeto list_df_matrix

tmp<-list_df_matrix[[1]]
matriz_1_to_add_noise<-as.matrix(tmp$normalized)
df0_no_noise<-as.data.frame(matriz_1_to_add_noise)
names(df0_no_noise)<-c("normalized")
matriz1_with_noise<-add.Gaussian.noise(matriz_1_to_add_noise, mean = 0, stddev = 0.01, symm = FALSE)
df1_with_noise<-as.data.frame(matriz1_with_noise)
names(df1_with_noise)<-c("normalized")

df1_with_noise<-cbind(tmp$M,df1_with_noise)
names(df1_with_noise)<-c("M","normalized")
df1_with_noise$normalized[which(df1_with_noise$normalized < 0)]<-0

#MAKE THE WAVELET WITH NOISE Nº 1 and identification

list_df_matrix<-NULL
list_df_matrix[[1]]<-df1_with_noise

cat("transform to wavelets", "\n")

wavelets_c<-cmpfun(wavelets)
wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$wavelet_mgf<-wavelets

df_spectrums_mgf$idx<-seq(1:(nrow(df_spectrums_mgf)))

cat("merge","\n")
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)

cat("calculating distances betweeen cosines", "\n")

dist_wav_c<-cmpfun(Dist_wav)
dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x), future.seed=TRUE)
df_results_tmp2$Dist<-dist_list
df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,Dist)),]
df_results_tmp2$wavelet<-NULL
df_results_tmp2$wavelet_mgf<-NULL

best_number<-10
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

cat("saving results", "\n")
new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")
new_name2<-paste(new_name,"dist",sep="_")
new_name3<-paste(new_name2,wave,sep="_")
new_name4<-paste0(new_name3,"NOISE_0.01_wp")

save(df_results_tmp2, file = paste0(new_name4,".Rdata"))

######### Calculate SNR to signal with noise nº 1

matriz2<-data.matrix(df1_with_noise$normalized)

SNR.Th <- 2
peakInfo <- peakDetectionCWT(matriz2,scales=1:30, SNR.Th=SNR.Th)
majorPeakInfo = peakInfo$majorPeakInfo
peakIndex <- majorPeakInfo$peakIndex

list1<-paste(peakInfo$majorPeakInfo$peakSNR)   #este identifica
distancia_id1<-as.numeric(df_results_tmp2$Dist[1])  #3.37

snr1<-data.frame("sum(peakSNR)"=sum(as.numeric(list1)),"dist_cosines"=distancia_id1, "signal"=1)

results_to_plot<-rbind(results_to_plot, snr1)

#################################################################################
#### añado más ruido a partir de la señal ruidosa 1

matriz2_with_noise<-add.Gaussian.noise(matriz_1_to_add_noise, mean = 0, stddev = 0.02, symm = FALSE)

df2_with_noise<-as.data.frame(matriz2_with_noise)
names(df2_with_noise)<-c("normalized")

df2_with_noise<-cbind(df1_with_noise$M,df2_with_noise)
names(df2_with_noise)<-c("M","normalized")
df2_with_noise$normalized[which(df2_with_noise$normalized < 0)]<-0

#MAKE THE WAVELET WITH NOISE Nº 1 and identification

list_df_matrix<-NULL
list_df_matrix[[1]]<-df2_with_noise

cat("transform to wavelets", "\n")

wavelets_c<-cmpfun(wavelets)
wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$wavelet_mgf<-wavelets

df_spectrums_mgf$idx<-seq(1:(nrow(df_spectrums_mgf)))

cat("merge","\n")
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)

cat("calculating distances betweeen cosines", "\n")

dist_wav_c<-cmpfun(Dist_wav)
dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x), future.seed=TRUE)
df_results_tmp2$Dist<-dist_list
df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,Dist)),]
df_results_tmp2$wavelet<-NULL
df_results_tmp2$wavelet_mgf<-NULL

best_number<-10
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

cat("saving results", "\n")
new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")
new_name2<-paste(new_name,"dist",sep="_")
new_name3<-paste(new_name2,wave,sep="_")
new_name4<-paste0(new_name3,"NOISE_0.02_wp")

save(df_results_tmp2, file = paste0(new_name4,".Rdata"))

######### Calculate SNR to signal with noise nº 2

matriz3<-data.matrix(df2_with_noise$normalized)

SNR.Th <- 2
peakInfo <- peakDetectionCWT(matriz3,scales=1:30, SNR.Th=SNR.Th)
majorPeakInfo = peakInfo$majorPeakInfo
peakIndex <- majorPeakInfo$peakIndex

list2<-paste(peakInfo$majorPeakInfo$peakSNR)   #este identifica

distancia_id2<-as.numeric(df_results_tmp2$Dist[1])

snr2<-data.frame("sum(peakSNR)"=sum(as.numeric(list2)),"dist_cosines"=distancia_id2, "signal"=2)

results_to_plot<-rbind(results_to_plot, snr2)

#################################################################################################
#################################################################################################

#################################################################################
#### añado más ruido a partir de la señal ruidosa 2

matriz3_with_noise<-add.Gaussian.noise(matriz_1_to_add_noise, mean = 0, stddev = 0.03, symm = FALSE)

df3_with_noise<-as.data.frame(matriz3_with_noise)
names(df3_with_noise)<-c("normalized")

df3_with_noise<-cbind(df1_with_noise$M,df3_with_noise)
names(df3_with_noise)<-c("M","normalized")
df3_with_noise$normalized[which(df3_with_noise$normalized < 0)]<-0

#MAKE THE WAVELET WITH NOISE Nº 3 and identification

list_df_matrix<-NULL
list_df_matrix[[1]]<-df3_with_noise

cat("transform to wavelets", "\n")

wavelets_c<-cmpfun(wavelets)
wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$wavelet_mgf<-wavelets

df_spectrums_mgf$idx<-seq(1:(nrow(df_spectrums_mgf)))

cat("merge","\n")
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)

cat("calculating distances betweeen cosines", "\n")

dist_wav_c<-cmpfun(Dist_wav)
dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x), future.seed=TRUE)
df_results_tmp2$Dist<-dist_list
df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,Dist)),]
df_results_tmp2$wavelet<-NULL
df_results_tmp2$wavelet_mgf<-NULL

best_number<-10
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

cat("saving results", "\n")
new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")
new_name2<-paste(new_name,"dist",sep="_")
new_name3<-paste(new_name2,wave,sep="_")
new_name4<-paste0(new_name3,"NOISE_0.03_wp")

save(df_results_tmp2, file = paste0(new_name4,".Rdata"))

######### Calculate SNR to signal with noise nº 3

matriz4<-data.matrix(df3_with_noise$normalized)

SNR.Th <- 2
peakInfo <- peakDetectionCWT(matriz4,scales=1:30, SNR.Th=SNR.Th)
majorPeakInfo = peakInfo$majorPeakInfo
peakIndex <- majorPeakInfo$peakIndex

list3<-paste(peakInfo$majorPeakInfo$peakSNR)   #este identifica

distancia_id3<-as.numeric(df_results_tmp2$Dist[1])

snr3<-data.frame("sum(peakSNR)"=sum(as.numeric(list3)),"dist_cosines"=distancia_id3, "signal"=3)

results_to_plot<-rbind(results_to_plot, snr3)
#results_to_plot<-results_to_plot[-c(4),]

#### add more noise to df3_with_noise

matriz4_with_noise<-add.Gaussian.noise(matriz_1_to_add_noise, mean = 0, stddev = 0.1, symm = FALSE)

df4_with_noise<-as.data.frame(matriz4_with_noise)
names(df4_with_noise)<-c("normalized")

df4_with_noise<-cbind(df1_with_noise$M,df4_with_noise)
names(df4_with_noise)<-c("M","normalized")
df4_with_noise$normalized[which(df4_with_noise$normalized < 0)]<-0


list_df_matrix<-NULL
list_df_matrix[[1]]<-df4_with_noise

cat("transform to wavelets", "\n")

wavelets_c<-cmpfun(wavelets)
wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$wavelet_mgf<-wavelets

df_spectrums_mgf$idx<-seq(1:(nrow(df_spectrums_mgf)))

cat("merge","\n")
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)

cat("calculating distances betweeen cosines", "\n")

dist_wav_c<-cmpfun(Dist_wav)
dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x), future.seed=TRUE)
df_results_tmp2$Dist<-dist_list
df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,Dist)),]
df_results_tmp2$wavelet<-NULL
df_results_tmp2$wavelet_mgf<-NULL

best_number<-10
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

cat("saving results", "\n")
new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")
new_name2<-paste(new_name,"dist",sep="_")
new_name3<-paste(new_name2,wave,sep="_")
new_name4<-paste0(new_name3,"NOISE_0.1_wp")

save(df_results_tmp2, file = paste0(new_name4,".Rdata"))

######### Calculate SNR to signal with noise nº 3

matriz5<-data.matrix(df4_with_noise$normalized)

SNR.Th <- 2
peakInfo <- peakDetectionCWT(matriz5,scales=1:30, SNR.Th=SNR.Th)
majorPeakInfo = peakInfo$majorPeakInfo
peakIndex <- majorPeakInfo$peakIndex

list4<-paste(peakInfo$majorPeakInfo$peakSNR)   #este identifica

distancia_id4<-as.numeric(df_results_tmp2$Dist[1])

snr4<-data.frame("sum(peakSNR)"=sum(as.numeric(list4)),"dist_cosines"=distancia_id4, "signal"=4)

library(ggplot2)

results_to_plot<-rbind(results_to_plot, snr4)
#results_to_plot<-results_to_plot[-c(4),]

##################################### add more noise to latest signal

matriz5_with_noise<-add.Gaussian.noise(matriz4_with_noise, mean = 0, stddev = 1, symm = FALSE)

df5_with_noise<-as.data.frame(matriz5_with_noise)
names(df5_with_noise)<-c("normalized")

df5_with_noise<-cbind(df1_with_noise$M,df5_with_noise)
names(df5_with_noise)<-c("M","normalized")
df5_with_noise$normalized[which(df5_with_noise$normalized < 0)]<-0


list_df_matrix<-NULL
list_df_matrix[[1]]<-df5_with_noise

cat("transform to wavelets", "\n")

wavelets_c<-cmpfun(wavelets)
wavelets<-future_lapply(list_df_matrix, FUN = function(x) wavelets_c(x), future.seed=TRUE)

df_spectrums_mgf$peaks<-NULL
df_spectrums_mgf$wavelet_mgf<-wavelets

df_spectrums_mgf$idx<-seq(1:(nrow(df_spectrums_mgf)))

cat("merge","\n")
df_results_tmp<-merge(df_spectrums_mgf,df_spectrums_selected_for_each_spectrum, by="idx", all.x=TRUE )

#eliminar no hits
df_results_tmp2<-df_results_tmp[which(!is.na(df_results_tmp$Name)),]
rm(df_results_tmp)

cat("calculating distances betweeen cosines", "\n")

dist_wav_c<-cmpfun(Dist_wav)
dist_list<-future_apply(df_results_tmp2,1, FUN = function(x) dist_wav_c(x), future.seed=TRUE)
df_results_tmp2$Dist<-dist_list
df_results_tmp2<-df_results_tmp2[with(df_results_tmp2, order(idx,Dist)),]
df_results_tmp2$wavelet<-NULL
df_results_tmp2$wavelet_mgf<-NULL

best_number<-10
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

cat("saving results", "\n")
new_name<-paste(substr(mgf, 1, nchar(mgf) - 3), "results", sep="")
new_name2<-paste(new_name,"dist",sep="_")
new_name3<-paste(new_name2,wave,sep="_")
new_name4<-paste0(new_name3,"NOISE_1_wp")

save(df_results_tmp2, file = paste0(new_name4,".Rdata"))

matriz6<-data.matrix(df5_with_noise$normalized)

SNR.Th <- 2
peakInfo <- peakDetectionCWT(matriz6,scales=1:30, SNR.Th=SNR.Th)
majorPeakInfo = peakInfo$majorPeakInfo
peakIndex <- majorPeakInfo$peakIndex

list5<-paste(peakInfo$majorPeakInfo$peakSNR)   #este identifica

distancia_id5<-as.numeric(df_results_tmp2$Dist[1])

snr5<-data.frame("sum(peakSNR)"=sum(as.numeric(list5)),"dist_cosines"=distancia_id5, "signal"=6)

results_to_plot<-rbind(results_to_plot, snr5)




###################### PLOTS ######################################################################

pdf("/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/plots_wavelet_SNR_DIST_Evolution.pdf")
ggplot()+
    geom_line(data=results_to_plot, aes(x=signal, y=sum.peakSNR.), color="red") +
    geom_line(data=results_to_plot, aes(x=signal, y=dist_cosines ), color="green")+
    labs(x="signals", y="SNR (red) & dist_cosines (green)")
dev.off()

pdf("/home/margaret/data/pepe/12_SIMULACIONES_WAVESPECTRO/señales_con_ruido.pdf")
par(mfrow=c(3,2))
plot(df0_no_noise[[1]]$normalized, type="l")
plot(df1_with_noise$normalized, type="l")
plot(df2_with_noise$normalized, type="l")
plot(df3_with_noise$normalized, type="l")
plot(df4_with_noise$normalized, type="l")
dev.off()
