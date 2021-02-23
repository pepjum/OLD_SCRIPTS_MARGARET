#!/usr/bin/Rscript

args=(commandArgs(TRUE))

database<-args[1]   #sptxt de spectraST
spread_value<-args[2] #0   valor para máqunas modernas es 0
Dalton<-args[3]   #1
n_peaks<-args[4]   #filtering spectra by number of peaks
intensity<-args[5] #minimun intensity threshold

source("/home/margaret/data/01_Rscripts/A_Funciones/funcionesShotgun.R")

#mgf_dir="."
options(stringsAsFactors = FALSE)

library(stringr)
library(stringi)
library(wavethresh)

old<-Sys.time()
print(old)

#spread_value<-0


readSpectraSTdb<-function(x){
    con <- file(x, open="r")
    db_raw <- readLines(con)
    close(con)

#db_dig_protsAll <- db_dig_raw[grep("Protein Name", db_dig_raw)]
    Names_spectra_ALL<-db_raw[grep("^Name",db_raw)]
    lib_id_spectra_ALL<-db_raw[grep("^LibID",db_raw)]
    MW_ALL<-db_raw[grep("^MW", db_raw)]
    PrecursorMZ_ALL<-db_raw[grep("^PrecursorMZ",db_raw)]
    Status_ALL<-db_raw[grep("^Status",db_raw)]
    Full_name_ALL<-db_raw[grep("^FullName",db_raw)]
    Comments_ALL<-db_raw[grep("^Comment",db_raw)]
    NumPeaks<-db_raw[grep("^NumPeaks",db_raw)]
    Names_spectra_ALL <- t(data.frame(strsplit(unlist(lapply(strsplit(Names_spectra_ALL, ": "), FUN = function(x) x[2])), "\\|")))

    colnames(Names_spectra_ALL) <- NULL
    rownames(Names_spectra_ALL) <- NULL
    colnames(Names_spectra_ALL)<-"Name"

    lib_id_spectra_ALL<-as.integer(str_extract(lib_id_spectra_ALL,'[0-9]+$'))
    PrecursorMZ_ALL<-as.numeric(str_extract(PrecursorMZ_ALL,'[0-9\\.]+'))
    Status_ALL<-sapply(strsplit(paste(Status_ALL)," "),'[',2)
    Full_name_ALL<-sapply(strsplit(paste(Full_name_ALL)," "),'[',2)
    NumPeaks<-as.numeric(str_extract(NumPeaks,'[0-9\\.]+'))
    MW_All<-as.numeric(str_extract(MW_ALL,'[0-9\\.]+'))


    tmp<-db_raw[-c(1:6)]  #elimino cabecera
    tmp<-lapply(tmp, FUN=function(x) x=paste0(x,"####")) #meto un fin de linea reconocible
    tmp<-lapply(tmp, FUN=function(x) if(grepl("^NumPeaks",x)){x=paste0(x,"@")}else{x=x})
    tmp<-stri_paste(tmp,collapse='') #uno todo en un string
    tmp<-strsplit(paste(tmp),"#Name")  #separo por name y asi separamos los espectros
    tmp <- unlist(lapply(tmp, FUN = function(x) x[x != ""]))
    tmp<-sapply(strsplit(paste(tmp),"@"),'[',3)

    SpectraST_df <- cbind(Names_spectra_ALL, data.frame("Lib_ID" = paste(lib_id_spectra_ALL), "MW" = paste(MW_All), "PrecursorMZ" = paste(PrecursorMZ_ALL), "Status" = paste(Status_ALL), "FullName" = paste(Full_name_ALL), "NumPeaks" = paste(NumPeaks), "Comments"=paste(Comments_ALL), "Peaks"=paste(tmp)))

    return(SpectraST_df)
}


#####CALL TO DB
cat("loading db","\n")
cat("Patience.... ","\n")
db_spectra<-readSpectraSTdb(paste(database))



filteringSpectraByNumberOfPeaks<-function(dataframe, n_of_peaks){
    selected_spectra_filtered<-dataframe[which(as.numeric(paste(dataframe$NumPeaks)) >= as.numeric(n_of_peaks)),]
    return(selected_spectra_filtered)
}
cat("filtering spectra by number of peaks","\n")
db_spectra_filtered<-filteringSpectraByNumberOfPeaks(db_spectra,n_peaks)


cat("filtering spectra by intensity","\n")
filteringSpectraByIntensity<-function(dataframe, min_mzvalue){
    selected_spectra_filtered_Mz<-dataframe[which(as.numeric(paste(dataframe$PrecursorMZ)) >= as.numeric(min_mzvalue)),]
    return(selected_spectra_filtered_Mz)
}

db_spectra_filtered_filtered_by_intensity<-filteringSpectraByIntensity(db_spectra_filtered, intensity)


list_all_peaks<-paste(db_spectra_filtered_filtered_by_intensity$Peaks)

db_spectra_filtered_filtered_by_intensity$Peaks<-NULL

bining<-function(dataframe,Dalton){
    dataframe<-as.data.frame(dataframe)
	new = 0
	old = nrow(dataframe)

	dataframe$M<-as.numeric(str_extract(paste(dataframe$M),'^[0-9]*'))
    dataframe$Z<-as.numeric(paste(dataframe$Z))
	while( old != new){

		tmp = data.frame('M' = NULL, 'Z' = NULL)
		SKIP_LINE = F
		iteration = nrow(dataframe)- 1
		for (i in 1:(iteration)){
			if(SKIP_LINE){
				SKIP_LINE = F
				next
			}
			dataframe<-dataframe[order(as.numeric(paste(dataframe$M))),]
			if(as.numeric(paste(dataframe[i+1, 'M'])) - as.numeric(paste(dataframe[i, 'M']))  <= as.numeric(Dalton)){
				#print('shrinking')
				maxM = max(as.numeric(paste(dataframe[i+1, 'M'])), as.numeric(paste(dataframe[i, 'M'] )))
				maxZ = max(as.numeric(paste(dataframe[i+1, 'Z'])), as.numeric(paste(dataframe[i, 'Z'] )))
				tmp = rbind(tmp, data.frame('M' = as.numeric(maxM), 'Z' = as.numeric(maxZ)))
				SKIP_LINE = T
			}else{
				tmp = rbind(tmp, dataframe[i,])
				SKIP_LINE = F
			}
		}
		tmp <- rbind(tmp, dataframe[i+1,])
		new <- nrow(tmp)
		old <- nrow(dataframe)
		dataframe <- tmp
	}
	if(as.numeric(paste(dataframe[nrow(dataframe), 'M'])) - as.numeric(paste(dataframe[nrow(dataframe)-1, 'M'])) <= as.numeric(Dalton)){
		maxM = max(as.numeric(paste(dataframe[nrow(dataframe), 'M'])), as.numeric(paste(dataframe[nrow(dataframe)-1 , 'M'] )))
		maxZ = max(as.numeric(paste(dataframe[nrow(dataframe), 'Z'])), as.numeric(paste(dataframe[nrow(dataframe)-1, 'Z'] )))
		dataframe2 <- dataframe[c(seq(1,nrow(dataframe)-2)), ]
		dataframe2 <- rbind(dataframe2, data.frame('M' = as.numeric(maxM), 'Z' = as.numeric(maxZ)))
	}else{
        dataframe2 <- dataframe

    }
    dataframe2<-unique(dataframe2)
	return(dataframe2)

}

####transform peaks info in dataframes to work in it
cat("transform peaks info in dataframes to work in it","\n")
list_df_peaks_not_binned<-list()
for(j in 1:length(list_all_peaks)){
    #cat(j, "of",length(list_all_peaks),"\n")
    tmp<-strsplit(paste(list_all_peaks[j]), "####")
    tmp<-tmp[!grepl('^###',tmp)]
    tmp<-lapply(tmp, FUN= function(x) x[!grepl('^###',x)])
    m<-paste(lapply(strsplit(unlist(tmp),"\\t"),'[',1))
    z<-paste(lapply(strsplit(unlist(tmp),"\\t"),'[',2))
    df_peaks<-data.frame("M"=m,"Z"=z)
    print(j)
    print(dim(df_peaks))
    list_df_peaks_not_binned[[j]]<-df_peaks

}
cat("saving database not binned","\n")
save(list_df_peaks_not_binned, file=paste0(database,"_NOT_BINNED.Rdata"))

cat("bining peaks","\n")
list_df_peaks_binned<-list()
for (k in 1: length(list_df_peaks_not_binned)){
    #cat(k,"\n")
    df_peaks_binned<-bining(list_df_peaks_not_binned[k],1)
    list_df_peaks_binned[[k]]<-df_peaks_binned
}
cat("saving database binned","\n")
save(list_df_peaks_binned,file=paste0(database,"_BINNED.Rdata"))

cat("spread out","\n")
list_df_peaks_binned_and_spread_out<-list()
for (i in 1:length(list_df_peaks_binned)){
    #cat(i, "of ", length(list_df_peaks_binned), "\n" )
    df<-as.data.frame(list_df_peaks_binned[i][[1]])
    df$percentage<-0

    df$percentage<-(as.numeric(paste(df$Z)))*as.numeric(paste(spread_value))/2    #el porcentaje del pico a transferir hay que dividirlo entre sus dos vecinos
    ###ahora quito el porcentaje de los picos originales
    df$peak_less<-as.numeric(paste(df$Z))-as.numeric(paste(df$percentage))    #el porcentaje del pico a transferir hay que dividirlo entre sus dos vecinos

    #sumar los pedazos de los adyacentes
    df$peaks_final<-0
    df$peaks_final[1]<-as.numeric(paste(df$peak_less[1]))+as.numeric(paste(df$percentage[2]))
    df$peaks_final[nrow(df)]<-as.numeric(paste(df$peak_less[nrow(df)]))+as.numeric(paste(df$percentage[nrow(df)-1]))

    for (j in 2:(nrow(df)-1)){
        df$peaks_final[j]<-as.numeric(paste(df$peak_less[j]))+as.numeric(paste(df$percentage[j-1]))+as.numeric(paste(df$percentage[j+1]))    #el porcentaje del pico a transferir hay que dividirlo entre sus dos vecinos
    }
    df_final<-df[,-c(2,3,4)]
    #df_final<-df
    #print(i)
    #print(dim(df_final))
    list_df_peaks_binned_and_spread_out[[i]]<-df_final
}
cat("saving database binned and 'spread out'","\n")
save(list_df_peaks_binned_and_spread_out, file=paste0(database,"_SPREAD_OUT.Rdata"))

#select range of spectrums before normalizing
list_selected<-list()
for (i in 1:length(list_df_peaks_binned_and_spread_out)){
    #cat(i, "of ", length(list_df_peaks_binned_and_spread_out), "\n" )
    df<-as.data.frame(list_df_peaks_binned_and_spread_out[i][[1]])
    df_ranges_tmp<-df[as.numeric(paste(df$M))>=500,]
    df_ranges_final<-df_ranges_tmp[as.numeric(paste(df_ranges_tmp$M))<=1500,]
    list_selected[[i]]<-df_ranges_final
}

cat("normalizing peaks","\n")
list_df_normalized<-list()
for (i in 1:length(list_selected)){
    #cat(i, "of ", length(list_selected), "\n" )
    df<-as.data.frame(list_selected[i][[1]])
    df$normalized<-as.numeric(df$peaks_final)/sqrt(sum(as.numeric(df$peaks_final)^2))
    df<-df[,c(1,3)]
    list_df_normalized[[i]]<-df
}
cat("saving database normalized only peaks","\n")
save(list_df_normalized, file=paste0(database,"_NORMALIZED.Rdata"))

list_df_matrix<-list()
for (i in 1:length(list_df_normalized)){
    #cat(i, "of ", length(list_df_normalized), "\n" )
    df<-as.data.frame(list_df_normalized[i][[1]])
    matrix<-data.frame(M=as.numeric(seq(1:1500)),normalized=0)
    if (nrow(df)==0){
        df<- matrix
    }
    matrix[which(as.numeric(matrix$M) %in% as.numeric(df$M)),]<-df$normalized
    matrix$M<-seq(1:nrow(matrix))
    df_ranges_tmp<-matrix[as.numeric(paste(matrix$M))>=500,]
    df_ranges_final<-df_ranges_tmp[as.numeric(paste(df_ranges_tmp$M))<=1500,]

    list_df_matrix[[i]]<-df_ranges_final
}



cat("transform dataframes into a list and rebuild  final database","\n")

list_linecod<-c()
for (i in 1:length(list_df_matrix)){
    #cat(i,"/",length(list_df_matrix),"\n")
    df<-as.data.frame(list_df_matrix[i][[1]])
    df$line<-paste0(df$M,"\t",df$normalized,"###")
    linecod<-paste0(df$line)
    linecod<-paste(linecod,collapse="")
    list_linecod<-c(list_linecod, linecod)
}

db_spectra_filtered_filtered_by_intensity$Peaks<-paste(list_linecod)

save(db_spectra_filtered_filtered_by_intensity, file=paste0(database,"_range500-1500.Rdata"))
write_feather(db_spectra_filtered_filtered_by_intensity, file=paste0(database,"_range500-1500.feather"))


### wavelets
wavelets<-list()
for (i in 1: length(list_df_matrix)){
    #cat (i, "/", length(list_df_matrix), "\n")
    df <- as.data.frame(list_df_matrix[i][[1]])
    red<-c(as.list(df$normalized),rep(0,23))
    wavelets[[i]]<- wd(as.numeric(paste(red)), filter.number=10, family ="DaubLeAsymm", type="wavelet")

}

# unir aproximacion y detalle en una misma lista
wavelets_vectors_A_D<-list()
for (i in 1: length(wavelets)){
    #cat (i, "/", length(wavelets), "\n")
    wavelets_vectors_A_D[[i]]<- c(wavelets[i][[1]]$C,wavelets[i][[1]]$C )

}


db_spectra_filtered_filtered_by_intensity$wavelet<-paste(wavelets_vectors_A_D)
save(db_spectra_filtered_filtered_by_intensity, file=paste0(database,"_range500-1500_wavelet.Rdata"))
write_feather(db_spectra_filtered_filtered_by_intensity, file=paste0(database,"_range500-1500_wavelet.feather"))
print(new)
cat("end of processing database")
