
selecting_spectrums<-function(df_spectrums_mgf,range,db_spectra){
        PrecursorMGF<-as.numeric(paste(df_spectrums_mgf["mz_value"]))
        i<-as.numeric(paste(df_spectrums_mgf["idx"]))
        last<-selectSpectraByMZrange(db_spectra,PrecursorMGF,range, i)
}

selecting_spectrums_v2<-function(df_spectrums_mgf,range,db_spectra){
        PrecursorMGF<-as.numeric(paste(df_spectrums_mgf["mz_value"]))
        i<-as.numeric(paste(df_spectrums_mgf["idx"]))
        last<-selectSpectraByMZrange_V2(db_spectra,PrecursorMGF,range, i)
}



selecting_spectrums_db<-function(df_spectrums_mgf,range,db_spectra){
        PrecursorMGF<-as.numeric(paste(df_spectrums_mgf["mz_value"]))
        i<-as.numeric(paste(df_spectrums_mgf["idx"]))
        last<-selectSpectraByMW(db_spectra,PrecursorMGF,range, i)
}


filteringSpectraByNumberOfPeaks<-function(dataframe, n_of_peaks){
    selected_spectra_filtered<-dataframe[which(as.numeric(paste(dataframe$NumPeaks)) >= as.numeric(n_of_peaks)),]
    return(selected_spectra_filtered)
}

filteringSpectraByIntensity<-function(dataframe, min_mzvalue){
    selected_spectra_filtered_Mz<-dataframe[which(as.numeric(paste(dataframe$MW)) >= as.numeric(min_mzvalue)),]
    return(selected_spectra_filtered_Mz)
}

selectSpectraByMZrange<-function(SpectraST_df,PrecursorMGF,range, ind){

    threshold_up <- (as.numeric(PrecursorMGF)+ as.numeric(range))
    threshold_down <- (as.numeric(PrecursorMGF)- as.numeric(range))
    l_up <- which(as.numeric(SpectraST_df$PrecursorMZ) < threshold_up)
    l_down <- which(as.numeric(SpectraST_df$PrecursorMZ) > threshold_down)
    tmp <- SpectraST_df[intersect(l_up,l_down ),]
    if(nrow(tmp)!=0){
        tmp$idx <- ind
        return(tmp)
    }else{
        return(tmp)
    }
}

selectSpectraByMZrange_V2<-function(SpectraST_df,PrecursorMGF,range, ind){

    threshold_up <- (as.numeric(PrecursorMGF)+ as.numeric(range))
    threshold_down <- (as.numeric(PrecursorMGF)- as.numeric(range))
    l_up <- which(as.numeric(SpectraST_df$PrecursorMZ) < threshold_up)
    l_down <- which(as.numeric(SpectraST_df$PrecursorMZ) > threshold_down)
    tmp <- SpectraST_df[intersect(l_up,l_down ),]
    tmp<-tmp[,c("Name","Lib_ID")]
    if(nrow(tmp)!=0){
        tmp$idx <- ind
        return(tmp)
    }else{
        return(tmp)
    }
}



selectSpectraByMW<-function(SpectraST_df,PrecursorMGF,range, ind){

    threshold_up <- (as.numeric(PrecursorMGF)+ as.numeric(range))
    threshold_down <- (as.numeric(PrecursorMGF)- as.numeric(range))
    l_up <- which(as.numeric(SpectraST_df$MW) < threshold_up)
    l_down <- which(as.numeric(SpectraST_df$MW) > threshold_down)
    tmp <- SpectraST_df[intersect(l_up,l_down ),]
    if(nrow(tmp)!=0){
        tmp$idx <- ind
        return(tmp)
    }else{
        return(tmp)
    }
}



bining<-function(dataframe,Dalton){
    dataframe<-as.data.frame(dataframe)
    dataframe<-dataframe[!(dataframe$M=="NA"),]
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


bining2<-function(dataframe_lista,Dalton){
    dataframe_list<-list()
    for (i in 1:length(dataframe_lista)){
        cat(i,"\n")
        dataframe<-as.data.frame(dataframe_lista[i])
        dataframe<-na.omit(dataframe)
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

        dataframe_list[[i]]<-dataframe2
    }
	return(dataframe_list)

}






bining_aprox<-function(dataframe){
    dataframe<-as.data.frame(dataframe)
	output<-approx(dataframe["M"],dataframe["Z"])
}


convert_df<-function(list_all_peaks_index){
    tmp<-strsplit(paste(list_all_peaks_index), "##")
    tmp<-paste(lapply(strsplit(unlist(tmp)," "),'[',2))
    m<-paste(lapply(strsplit(unlist(tmp),"\t"),'[',1))
    z<-paste(lapply(strsplit(unlist(tmp),"\t"),'[',2))
    df_peaks<-data.frame("M"=m,"Z"=z)
}

convert_df3<-function(list_all_peaks_index){
    tmp<-strsplit(paste(list_all_peaks_index), "###")
    m<-paste(lapply(strsplit(unlist(tmp)," "),'[',1))
    z<-paste(lapply(strsplit(unlist(tmp)," "),'[',2))
    df_peaks<-data.frame("M"=m,"Z"=z)
}

convert_df_2<-function(list_all_peaks_index){
    tmp<-strsplit(paste(list_all_peaks_index), "####")
    tmp<-tmp[!grepl('^###',tmp)]
    tmp<-lapply(tmp, FUN= function(x) x[!grepl('^###',x)])
    m<-paste(lapply(strsplit(unlist(tmp),"\\t"),'[',1))
    z<-paste(lapply(strsplit(unlist(tmp),"\\t"),'[',2))
    df_peaks<-data.frame("M"=m,"Z"=z)
}


spread_out<-function(list,value){
    list_df_peaks_binned_and_spread_out<-list()
    for (i in 1:length(list)){
        if (i %% 100 ==0){
            #cat(i, "of ", length(list_df_peaks_binned), "\n" )
        }
        df<-as.data.frame(list[i][[1]])
        df$percentage<-0

        df$percentage<-(as.numeric(paste(df$Z)))*as.numeric(paste(value))/2    #el porcentaje del pico a transferir hay que dividirlo entre sus dos vecinos
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
        #print(dim(df_final))
        names(df_final)[2]<-"Z"
        list_df_peaks_binned_and_spread_out[[i]]<-df_final
        return(list_df_peaks_binned_and_spread_out)
    }
}

select_range<-function(list){
    df<-as.data.frame(list)
    df_ranges_tmp<-df[as.numeric(paste(df$M))>=150,]
    df_ranges_final<-df_ranges_tmp[as.numeric(paste(df_ranges_tmp$M))<=1500,]
}

normalization<-function(list){
    df<-as.data.frame(list)
    df$normalized<-as.numeric(df$Z)/sqrt(sum(as.numeric(df$Z)^2))
    df<-df[,c(1,3)]
}

normalization_w<-function(list){
    wavelet<-unlist(list)
    values_peaks<-wavelet[(wavelet !=0)]
    positions_peaks<-which(wavelet !=0)
    matrix<-data.frame("m"=seq(1:1024),"n"=0)
    w_norm<- values_peaks - min(values_peaks)/(max(values_peaks)- min(values_peaks))
    #w_norm<-values_peaks/sqrt(sum(values_peaks^2))
    datos<-data.frame("m"=positions_peaks,"n"=w_norm)
    matrix[which(as.numeric(matrix$m) %in% as.numeric(datos$m)),]<-datos$n
    output<-matrix$n
}

normalization_w2<-function(list){
    wavelet<-unlist(list)
    w_norm<- wavelet/sqrt(sum(wavelet^2))
}



create_full_vector<-function(list){
    df<-as.data.frame(list)
    matrix<-data.frame(M=as.numeric(seq(1:1500)),normalized=0)
    if (nrow(df)==0){
        df<- matrix
    }
    matrix[which(as.numeric(matrix$M) %in% as.numeric(df$M)),]<-df$normalized
    matrix$M<-seq(1:nrow(matrix))
    df_ranges_tmp<-matrix[as.numeric(paste(matrix$M))>=150,]
    df_ranges_final<-df_ranges_tmp[as.numeric(paste(df_ranges_tmp$M))<=1500,]

}

codification<-function(list){
    df<-as.data.frame(list)
    df$line<-paste0(df$M,"\t",df$normalized,"###")
    linecod<-paste0(df$line)
    linecod<-paste(linecod,collapse="")
}

Dot_and_DB<-function(df_results_tmp2){
    tmp_mgf<-strsplit(paste(df_results_tmp2["Peaks_mgf"]), "###")
    m_mgf<-paste(lapply(strsplit(unlist(tmp_mgf),"\t"),'[',1))
    z_mgf<-paste(lapply(strsplit(unlist(tmp_mgf),"\t"),'[',2))
    df_peaks_mgf<-data.frame("M"=as.numeric(m_mgf),"Z"=as.numeric(z_mgf))
    tmp_database<-strsplit(paste(df_results_tmp2["Peaks"]), "###")
    m_database<-paste(lapply(strsplit(unlist(tmp_database),"\t"),'[',1))
    z_database<-paste(lapply(strsplit(unlist(tmp_database),"\t"),'[',2))
    df_peaks_database<-data.frame("M"=as.numeric(m_database),"Z"=as.numeric(z_database))
    df_peaks_mgf_tmp1<-df_peaks_mgf[(df_peaks_mgf$M >150),]
    df_peaks_mgf_tmp2<-df_peaks_mgf_tmp1[(df_peaks_mgf_tmp1$M <1500),]
    df_peaks_database_tmp1<-df_peaks_database[(df_peaks_database$M >150),]
    df_peaks_database_tmp2<-df_peaks_database_tmp1[(df_peaks_database_tmp1$M <1500),]
    D<-sum(as.numeric(paste(df_peaks_mgf_tmp2$Z))*(as.numeric(paste(df_peaks_database_tmp2$Z))))
    mgf_square<-as.numeric(df_peaks_mgf_tmp2$Z)^2
    peaks_square<-as.numeric(df_peaks_database_tmp2$Z)^2
    sum_dot<-sum((as.numeric(mgf_square)) * (as.numeric(peaks_square)))
    DB<-sqrt(as.numeric(sum_dot))/as.numeric(D)
    res<-data.frame("D"= D, "DB"=DB)
}

Dot_and_DB_2<-function(df_results_tmp2, db_spectra){
    tmp_mgf<-strsplit(paste(df_results_tmp2["Peaks_mgf"]), "###")
    m_mgf<-paste(lapply(strsplit(unlist(tmp_mgf),"\t"),'[',1))
    z_mgf<-paste(lapply(strsplit(unlist(tmp_mgf),"\t"),'[',2))
    df_peaks_mgf<-data.frame("M"=as.numeric(m_mgf),"Z"=as.numeric(z_mgf))

    Peaks_database<-db_spectra[which(db_spectra$Lib_ID == df_results_tmp2$Lib_ID), "Peaks"]
    tmp_database<-strsplit(paste(Peaks_database), "###")
    m_database<-paste(lapply(strsplit(unlist(tmp_database),"\t"),'[',1))
    z_database<-paste(lapply(strsplit(unlist(tmp_database),"\t"),'[',2))
    df_peaks_database<-data.frame("M"=as.numeric(m_database),"Z"=as.numeric(z_database))
    df_peaks_mgf_tmp1<-df_peaks_mgf[(df_peaks_mgf$M >150),]
    df_peaks_mgf_tmp2<-df_peaks_mgf_tmp1[(df_peaks_mgf_tmp1$M <1500),]
    df_peaks_database_tmp1<-df_peaks_database[(df_peaks_database$M >150),]
    df_peaks_database_tmp2<-df_peaks_database_tmp1[(df_peaks_database_tmp1$M <1500),]
    D<-sum(as.numeric(paste(df_peaks_mgf_tmp2$Z))*(as.numeric(paste(df_peaks_database_tmp2$Z))))
    mgf_square<-as.numeric(df_peaks_mgf_tmp2$Z)^2
    peaks_square<-as.numeric(df_peaks_database_tmp2$Z)^2
    sum_dot<-sum((as.numeric(mgf_square)) * (as.numeric(peaks_square)))
    DB<-sqrt(as.numeric(sum_dot))/as.numeric(D)
    res<-data.frame("D"= D, "DB"=DB)
}


Dot_and_DB_waw<-function(df_results_tmp2){
    wave_mgf<-unlist(df_results_tmp2["wavelet_mgf"])
    wave_db<-unlist(df_results_tmp2["wavelet"])
    D<-sum((wave_mgf)*(wave_db))
    mgf_square<-(wave_mgf)^2
    peaks_square<-(wave_db)^2
    sum_dot<-sum((mgf_square) * (peaks_square))
    DB<-sqrt(as.numeric(sum_dot))/as.numeric(D)
    res<-data.frame("D"= D, "DB"=DB)
}

Dot_and_DB_waw_2<-function(df_results_tmp2, db_spectra){
    wave_mgf<-unlist(df_results_tmp2["wavelet_mgf"])
    lib_id<-paste(df_results_tmp2["Lib_ID"])
    wave_db<-unlist(db_spectra[which(paste(db_spectra$Lib_ID) == lib_id),"wavelet"])

    D<-sum((wave_mgf)*(wave_db))
    mgf_square<-(wave_mgf)^2
    peaks_square<-(wave_db)^2
    sum_dot<-sum((mgf_square) * (peaks_square))
    DB<-sqrt(as.numeric(sum_dot))/as.numeric(D)
    res<-data.frame("D"= D, "DB"=DB)
}


Dist_wav<-function(df_results_tmp2){
    wave_mgf<-unlist(df_results_tmp2["wavelet_mgf"])
    wave_db<-unlist(df_results_tmp2["wavelet"])
    Dist<-dist(rbind(wave_mgf,wave_db),method="euclidean")
}

Dist_wav_v2<-function(df_results_tmp2, db_spectra){
    wave_mgf<-unlist(df_results_tmp2["wavelet_mgf"])
    lib_id<-paste(df_results_tmp2["Lib_ID"])
    wave_db<-unlist(db_spectra[which(paste(db_spectra$Lib_ID) == lib_id),"wavelet"])

    Dist<-dist(rbind(wave_mgf,wave_db),method="euclidean")
}


Dist_cosine_wav<-function(df_results_tmp2){
    library(stylo)
    wave_mgf<-unlist(df_results_tmp2["wavelet_mgf"])
    wave_db<-unlist(df_results_tmp2["wavelet"])
    dataset<-rbind(wave_mgf, wave_db)
    #colnames(dataset)<-c("spec","db_spec")
    dist_cosines<-dist.cosine(dataset)

}

Dist_cosine_wav_v2<-function(df_results_tmp2, db_spectra){
    library(stylo)
    wave_mgf<-unlist(df_results_tmp2["wavelet_mgf"])
    lib_id<-paste(df_results_tmp2["Lib_ID"])
    wave_db<-unlist(db_spectra[which(paste(db_spectra$Lib_ID) == lib_id),"wavelet"])
    dataset<-rbind(wave_mgf, wave_db)
    #colnames(dataset)<-c("spec","db_spec")
    dist_cosines<-dist.cosine(dataset)

}

Dist_cosine_wav_v3<-function(df_results_tmp2){
    library(proxy)
        wave_mgf<-unlist(df_results_tmp2["wavelet_mgf"])
    wave_db<-unlist(df_results_tmp2["wavelet"])
    dataset<-rbind(wave_mgf, wave_db)
    dist_cosines<-dist(dataset,method="cosine")
}

mutual_information<-function(df_results_tmp2){
    wave_mgf<-as.numeric(paste(unlist(df_results_tmp2["wavelet_mgf"][1])))
    wave_db<-as.numeric(paste(unlist(df_results_tmp2["wavelet"][1])))
    results<-JSD(wave_mgf,wave_db)
}

mutual_information_v2<-function(df_results_tmp2, db_spectra){
    wave_mgf<-unlist(df_results_tmp2["wavelet_mgf"])
    lib_id<-paste(df_results_tmp2["Lib_ID"])
    wave_db<-unlist(db_spectra[which(paste(db_spectra$Lib_ID) == lib_id),"wavelet"])

    results<-JSD(wave_mgf,wave_db)
}



JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(wave_mgf,wave_db) {
        tmp <- (wave_mgf/wave_db)
        xx <- tmp[tmp>0]
        xx<-xx[is.finite(xx)]
        sum(xx * log(xx))
}



F_score<-function(dataframe){
    if (as.numeric(dataframe["DB"]) < 0.1){
        dataframe["Fscore"]<- (0.6*(as.numeric(dataframe["D"])))+(0.4*(as.numeric(dataframe["Delta"]))) - 0.12
    }else if (as.numeric(dataframe["DB"]) > 0.35 && as.numeric(dataframe["DB"]) <= 0.4 ){
        dataframe["Fscore"]<- (0.6*(as.numeric(dataframe["D"])))+(0.4*(as.numeric(dataframe["Delta"]))) - 0.12
    }else if(as.numeric(dataframe["DB"]) > 0.4 && as.numeric(dataframe["DB"]) <= 0.45){
        dataframe["Fscore"]<- (0.6*(as.numeric(dataframe["D"])))+(0.4*(as.numeric(dataframe["Delta"]))) - 0.18
    }else if(as.numeric(dataframe["DB"]) > 0.45){
        dataframe["Fscore"]<- (0.6*(as.numeric(dataframe["D"])))+(0.4*(as.numeric(dataframe["Delta"]))) - 0.18
    }else{
        dataframe["Fscore"]<- (0.6*(as.numeric(dataframe["D"])))+(0.4*(as.numeric(dataframe["Delta"]))) - 0
    }
}


Delta<-function(dataframe){

    dataframe$Delta <- NA
    for (i in unique(dataframe$idx)){
        tmp<-dataframe[which(dataframe$idx == i),]
        sortedD <- sort(as.numeric(tmp$D), decreasing = T)
        if (nrow(tmp) >= 2){
            delta <- sortedD[1] - sortedD[2]
        }else{
            delta <- sortedD[1]
        }
        dataframe[which(dataframe$idx == i),'Delta'] <- delta
    }
    return(dataframe)
}

Delta2<-function(dataframe){
    list_dataframes<-list()
    for (i in unique(dataframe$idx)){
        tmp<-dataframe[which(dataframe$idx == i),]
        tmp_max1<-tmp[(as.numeric(paste(tmp$D)) ==max(as.numeric(paste(tmp$D)))),]
        tmp_max2<-tmp[-c(as.numeric(paste(tmp$D)) ==max(as.numeric(paste(tmp$D)))),]
        if (nrow(tmp_max2)==0){
            D3<-as.numeric(tmp_max1$D)
        }else{
            tmp_max3<-tmp_max2[(as.numeric(paste(tmp_max2$D)) ==max(as.numeric(paste(tmp_max2$D)))),]
            D1<-as.numeric(tmp_max1$D)
            D2<-as.numeric(tmp_max3$D)
            D3<-(D1-D2)
        }
        tmp[tmp$idx == i, 'Delta']<-as.numeric(D3)
        list_dataframes[[i]]<-tmp
    }
    return(list_dataframes)
}





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
    Names_spectra_ALL<-sapply(strsplit(paste(Names_spectra_ALL),": "),"[",2)

    lib_id_spectra_ALL<-as.integer(str_extract(lib_id_spectra_ALL,'[0-9]+$'))
    PrecursorMZ_ALL<-as.numeric(str_extract(PrecursorMZ_ALL,'[0-9\\.]+'))
    Status_ALL<-sapply(strsplit(paste(Status_ALL)," "),'[',2)
    Full_name_ALL<-sapply(strsplit(paste(Full_name_ALL)," "),'[',2)
    NumPeaks<-as.numeric(str_extract(NumPeaks,'[0-9\\.]+'))
    MW_All<-as.numeric(str_extract(MW_ALL,'[0-9\\.]+'))


    tmp<-db_raw[-c(1:6)]  #elimino cabecera
    tmp<-future_lapply(tmp, FUN=function(x) x=paste0(x,"####"),future.seed=TRUE) #meto un fin de linea reconocible
    tmp<-lapply(tmp, FUN=function(x) if(grepl("^NumPeaks",x)){x=paste0(x,"@")}else{x=x})
    tmp<-stri_paste(tmp,collapse='') #uno todo en un string
    tmp<-strsplit(paste(tmp),"#Name")  #separo por name y asi separamos los espectros
    tmp <- unlist(lapply(tmp, FUN = function(x) x[x != ""]))
    tmp<-sapply(strsplit(paste(tmp),"@"),'[',2)

    SpectraST_df <-data.frame("Name" = paste(Names_spectra_ALL), "Lib_ID" = paste(lib_id_spectra_ALL), "MW" = paste(MW_All), "PrecursorMZ" = paste(PrecursorMZ_ALL), "Status" = paste(Status_ALL), "FullName" = paste(Full_name_ALL), "NumPeaks" = paste(NumPeaks), "Comments"=paste(Comments_ALL), "Peaks"=paste(tmp))

    return(SpectraST_df)
}


readMSPFile<-function(x){
    con <- file(x, open="r")
    Names_spectra_ALL<-list()
    MW_ALL<-list()
    Comments_ALL<-list()
    NumPeaks_ALL<-list()
    Peaks_ALL<-list()
    PrecursorMZ_ALL<-list()
    counter<-0
    line <- 'Hello'
    reading <- T
    while (reading == T){
        line <- readLines(con, n=1)
        #print(reading)
        if(length(line) == 0){break}
        if(grepl("^Name", line)){
            counter = counter+1
            Names_spectra_ALL[counter]<-line
        }else if(grepl("^MW",line)){
            MW_ALL[counter]<-line
        }else if(grepl("^PrecursorMZ",line)){
            PrecursorMZ_ALL[counter]<-line
        }else if(grepl("^Comment",line)){
            Comments_ALL[counter]<-line
        }else if(grepl("^Num", line)){
            NumPeaks_ALL[counter]<-line
        }else if(grepl('^[0-9]', line)){
            if(is.null(Peaks_ALL[counter][1])){
                Peaks_ALL[counter]<-line
            }else{
                Peaks_ALL[counter]<-paste(Peaks_ALL[counter], '##',line)
            }
         }
    }
    close(con)
    Names_spectra_ALL<-sapply(strsplit(paste(Names_spectra_ALL),": "),"[",2)
    NumPeaks_ALL<-as.numeric(str_extract(NumPeaks_ALL,'[0-9\\.]+'))
    MW_All<-as.numeric(str_extract(MW_ALL,'[0-9\\.]+'))
    #Charge_all<-sapply(strsplit(paste(Names_spectra_ALL),"\\"),"[",2)
    #Charge_all<-as.numeric(str_extract(Charge_all,'[0-9\\.]+' ))
    PrecursorMZ_ALL<-as.numeric(str_extract(PrecursorMZ_ALL,'[0-9\\.]+'))
    Peaks_ALL<-sapply(strsplit(paste(Peaks_ALL),"NULL ##"),"[",2)

    ProteinNameTMP<-sapply(strsplit(paste(Comments_ALL),"NISTProtein="),"[",2)
    ProteinName<-sapply(strsplit(paste(ProteinNameTMP),"\\("),"[",1)
    #limpiar campos y generar df


    MSP_df <-data.frame("Name" = paste(Names_spectra_ALL), "Lib_ID" = paste(seq(1:length(Names_spectra_ALL))), "MW" = paste(MW_All), "PrecursorMZ" = paste(PrecursorMZ_ALL), "NumPeaks" = paste(NumPeaks_ALL), "Comments"=paste(Comments_ALL), "Peaks"=paste(Peaks_ALL), "ProteinName"=paste(ProteinName))
    return(MSP_df)

}







mgf_to_df<-function(mgf_file){
    name_spectrums<-c()
    mass_values<-c()
    peaks_list<-c()
    for (i in 1:length(mgf_file)){
        name_file<-mgf_file[[i]][2]
        name_spectrums<-c(name_spectrums,name_file)
        mass<-mgf_file[[i]][4]
        mass_values<-c(mass_values,mass)
        peaks<-mgf_file[[i]][6:length(mgf_file[[i]])-1]
        peaks_cod<-paste(peaks,collapse="###")
        peaks_list<-c(peaks_list, peaks_cod)

    }
    mass_values<-future_sapply(strsplit(paste(mass_values),"="),`[`,2, future.seed=TRUE)
    mass_values<-future_sapply(strsplit(paste(mass_values)," "),`[`,1, future.seed=TRUE)

    name_files<-future_sapply(strsplit(paste(name_spectrums),"="),`[`,2, future.seed=TRUE)
    name_spectrums<-future_sapply(strsplit(paste(name_files)," "),`[`,1, future.seed=TRUE)

    #creating df from mgf file

    df_spectrums_mgf<-data.frame("spectrum"=name_spectrums,"mz_value"=mass_values, "peaks"=paste(peaks_list))
    df_spectrums_mgf$idx<-seq(1:nrow(df_spectrums_mgf))
    return(df_spectrums_mgf)
}

mgf_to_df2<-function(mgf_file){
    name_spectrums<-c()
    mass_values<-c()
    peaks_list<-c()
    for (i in 1:length(mgf_file)){
        name_file<-mgf_file[[i]][5]
        name_spectrums<-c(name_spectrums,name_file)
        mass<-mgf_file[[i]][3]
        mass_values<-c(mass_values,mass)
        peaks<-mgf_file[[i]][7:length(mgf_file[[i]])-1]
        peaks<-gsub("\t"," ",peaks)
        peaks_cod<-paste(peaks,collapse="###")
        peaks_list<-c(peaks_list, peaks_cod)

    }
    mass_values<-future_sapply(strsplit(paste(mass_values),"="),`[`,2, future.seed=TRUE)
    mass_values<-future_sapply(strsplit(paste(mass_values)," "),`[`,1, future.seed=TRUE)
    #mass_values<-as.numeric(mass_values)+ 20

    name_files<-future_sapply(strsplit(paste(name_spectrums),"="),`[`,2, future.seed=TRUE)
    name_spectrums<-future_sapply(strsplit(paste(name_files)," "),`[`,1, future.seed=TRUE)

    df_spectrums_mgf<-data.frame("spectrum"=name_spectrums,"mz_value"=mass_values, "peaks"=paste(peaks_list))
    df_spectrums_mgf$idx<-seq(1:nrow(df_spectrums_mgf))
    return(df_spectrums_mgf)

}


mgf_to_df3<-function(mgf_file){
    name_spectrums<-c()
    mass_values<-c()
    peaks_list<-c()
    for (i in 1:length(mgf_file)){
        name_file<-mgf_file[[i]][2]
        name_spectrums<-c(name_spectrums,name_file)
        mass<-mgf_file[[i]][4]
        mass_values<-c(mass_values,mass)
        peaks<-mgf_file[[i]][6:length(mgf_file[[i]])]
        peaks_cod<-paste(peaks,collapse="###")
        peaks_list<-c(peaks_list, peaks_cod)

    }
    mass_values<-future_sapply(strsplit(paste(mass_values),"="),`[`,2, future.seed=TRUE)
    mass_values<-future_sapply(strsplit(paste(mass_values)," "),`[`,1, future.seed=TRUE)

    name_files<-future_sapply(strsplit(paste(name_spectrums),"="),`[`,2, future.seed=TRUE)
    name_spectrums<-future_sapply(strsplit(paste(name_files)," "),`[`,1, future.seed=TRUE)

    #creating df from mgf file

    df_spectrums_mgf<-data.frame("spectrum"=name_spectrums,"mz_value"=mass_values, "peaks"=paste(peaks_list))
    df_spectrums_mgf$idx<-seq(1:nrow(df_spectrums_mgf))
    return(df_spectrums_mgf)
}


mgf_to_df4<-function(mgf_file){
    name_spectrums<-c()
    mass_values<-c()
    peaks_list<-c()
    charge<-c()
    for (i in 1:length(mgf_file)){
        name_file<-mgf_file[[i]][2]
        name_spectrums<-c(name_spectrums,name_file)
        mass<-mgf_file[[i]][3]
        mass_values<-c(mass_values,mass)
        charges<-mgf_file[[i]][4]
        charge<-c(charge, charges)
        peaks<-mgf_file[[i]][5:(length(mgf_file[[i]])-1)]
        peaks_cod<-paste(peaks,collapse="###")
        peaks_list<-c(peaks_list, peaks_cod)

    }
    mass_values<-future_sapply(strsplit(paste(mass_values),"="),`[`,2, future.seed=TRUE)
    #mass_values<-future_sapply(strsplit(paste(mass_values)," "),`[`,1, future.seed=TRUE)

    name_files<-future_sapply(strsplit(paste(name_spectrums),"="),`[`,2, future.seed=TRUE)
    #name_spectrums<-future_sapply(strsplit(paste(name_files)," "),`[`,1, future.seed=TRUE)
    charges_tmp<-future_sapply(strsplit(paste(charge),"="),`[`,2, future.seed=TRUE)
    charges<-future_sapply(strsplit(paste(charges_tmp),"+"),`[`,1, future.seed=TRUE)
    #creating df from mgf file

    df_spectrums_mgf<-data.frame("spectrum"=name_files,"mz_value"=mass_values, "peaks"=paste(peaks_list), "charge"=paste(charges))
    df_spectrums_mgf$idx<-seq(1:nrow(df_spectrums_mgf))
    return(df_spectrums_mgf)
}






wavelets_for<-function(list){
    df<- as.data.frame(list)
    red<-c(as.list(df$normalized),rep(0,(2048-length(df$normalized))))   #zeropadding al final del vector para que de el minimo de tamaño permitido por la wavelet
    wavelet<-wd(as.numeric(paste(red)), filter.number=8, family ="DaubLeAsymm",bc="interval", verbose=FALSE, min.scale=1)
    wavelet_denoising<-threshold(wavelet, levels=1, policy="universal", by.level=FALSE, return.threshold=FALSE)
    wavelet_reconst<-wr(wavelet_denoising)
    wavelet_final<-wd(as.numeric(paste(wavelet_reconst)), filter.number=8, family ="DaubLeAsymm",bc="interval", verbose=FALSE, min.scale=1)
    wavelet_vector<-wavelet$transformed.vector
    #wavelet_sd<-sd(wavelet_vector)
}

wavelets_d<-function(list){
    df<- as.data.frame(list)
    red<-c(as.list(df$normalized),rep(0,(2048-length(df$normalized))))   #zeropadding al final del vector para que de el minimo de tamaño permitido por la wavelet
    wavelet<-wd(as.numeric(paste(red)), filter.number=10, family ="DaubLeAsymm")
    wavelet_vector<-c(unlist(wavelet["C"]),unlist(wavelet["D"]))
}


wavelets_packets<-function(list){
    df<- as.data.frame(list)
    red<-c(as.list(df$normalized),rep(0,(2048-length(df$normalized))))   #zeropadding al final del vector para que de el minimo de tamaño permitido por la wavelet
    wavelet<-wp(as.numeric(paste(red)), filter.number=10, family ="DaubLeAsymm", verbose=FALSE)
    wavesel<-wavelet$wp
    wavelet_vector<-unlist(split(wavesel, rep(1:ncol(wavesel), each = nrow(wavesel))))
    #wavelet_vector<-c(unlist(wavelet["C"]),unlist(wavelet["D"]))
}

wavelets_2<-function(list){
    df<- as.data.frame(list)
    red<-c(as.list(df$M),rep(0,(2048-length(df$Z))))   #zeropadding al final del vector para que de el minimo de tamaño permitido por la wavelet
    wavelet<-wp(as.numeric(paste(red)), filter.number=10, family ="DaubLeAsymm", verbose=FALSE)
    wavesel<-wavelet$wp
    wavelet_vector<-unlist(split(wavesel, rep(1:ncol(wavesel), each = nrow(wavesel))))
    #wavelet_vector<-c(unlist(wavelet["C"]),unlist(wavelet["D"]))
}


ReadMGFFile <- function(fileName)
{
	conn=file(fileName,open="r")
	mgfFile=readLines(conn)
	indBeginIons <- grep("BEGIN IONS", mgfFile)
	indEndIons <- grep("END IONS", mgfFile)
	indTitle <- indBeginIons + 1
	titlesMGF <- mgfFile[indTitle]
	titlesMGF <- strsplit(titlesMGF, "=")
	titlesMGF <- unlist(future_lapply(titlesMGF, FUN = function(x) x[2], future.seed=TRUE))
	mgfList <- future_lapply(1:length(indBeginIons), FUN = function(x) mgfFile[indBeginIons[x]:indEndIons[x]],future.seed=TRUE)
	names(mgfList) <- titlesMGF
	close(conn)
	return(mgfList)
}

ReadMGFFileMASSSimulator <- function(fileName)
{
	conn=file(fileName,open="r")
	mgfFile=readLines(conn)
	indBeginIons <- grep("BEGIN IONS", mgfFile)
	indEndIons <- grep("END IONS", mgfFile)
	indTitle <- indBeginIons + 1
    indPepMass <- indBeginIons + 2
    #indCharge <- indBeginIons + 3
	titlesMGF <- mgfFile[indTitle]
    PepMassMGF<- mgfFile[indPepMass]
	titlesMGF <- unlist(lapply(strsplit(titlesMGF, "="),"[",2))
	mgfList <- future_lapply(1:length(indBeginIons), FUN = function(x) mgfFile[indBeginIons[x]:indEndIons[x]],future.seed=TRUE)
	names(mgfList) <- titlesMGF
	close(conn)
	return(mgfList)
}
