#/usr/bin/Rscript
args=(commandArgs(TRUE))

directory<-args[1]
xmlFile_in <- args[2]

cat(system("date"),"\n")

#cat(xmlFile_in,"\n")

require(XML, quietly=T)
require(stringr,quietly=T)
require(gdata,quietly=T)
require(RCurl, quietly=T)
xmlFile_in<-lapply(strsplit(paste(xmlFile_in),"\\/"),"[",2)
#save(xml_data_raw, file = paste(gsub(".xml","",xmlFile), "_raw.rda", sep = ""))

parsedataXML<-function(xml_data_unlist){

    
    indexes <- which(names(xml_data_unlist) == "entry-list.entry.overview.protein-existence.value")
    cat("proteins =", length(indexes),"\n")

    xml_data_df <- data.frame("NX_code"=as.character(),
    "GO_code"=as.character()
    )

        for (i in 1:length(indexes)) {
        	if (i == 1) {
          # cat(i, "of ",length(indexes)-1,"\n")    
			index <- 1
			index_fin <- indexes[i] - 1
           # cat(index,"/",index_fin,"\n")
		    }else{
            #cat(i, "of ",length(indexes)-1,"\n")
			index <- indexes[i-1]
			index_fin <- indexes[i] - 1
            #cat(index,"/",index_fin,"\n")
		    }
            tmp <- xml_data_unlist[index:index_fin]
            tmp_df <- data.frame("name" = names(tmp), "value" = tmp)

            if(length(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"]))) ==1 & (length(unique(paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"]))) >=1)){
                
                subset_GO_tmp<-paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"])
                subset_KW<-subset_GO_tmp[startsWith(subset_GO_tmp, "KW")]
                if(length(subset_KW)!=0){
                    output_KW<-c()
                   for(KW in subset_KW){
                       #cat("KW ", KW,"\n"
                       sys.sleep(2)
                       file<-readLines(paste0("https://www.uniprot.org/keywords/",KW), warn=F)
                       GO_out<-str_extract(grep("QuickGO", file, value=T),'GO:[0-9]+')
                        if(length(GO_out)==1){
                            output_KW<-c(output_KW, GO_out)
                        }else if(length(GO_out) > 1){
                            #cat("KW distinto de 1",KW, "\n")
                        }else{
                            next
                        } 
                     
                   }

                }
                subset_GO<-subset_GO_tmp[startsWith(subset_GO_tmp, "GO:")]
                subset_GO<-unique(subset_GO)
                if(exists("output_KW")){
                    if(length(subset_GO)==0 & length(output_KW)==0){
                        output_total<-"NA"
                    }else if(length(subset_GO) >0 & length(output_KW) >0){
                        output_total<-c(subset_GO, output_KW)
                        output_total<-unique(output_total)
                    }else if(length(subset_GO) >0 & length(output_KW)==0){
                        output_total<-subset_GO
                    }else if(length(subset_GO)==0 & length(output_KW)>0){
                        output_total<-unique(output_KW)
                    } 
                }else{
                    if(length(subset_GO)==0){
                        output_total<-"NA"
                    }else{
                    output_total<-subset_GO
                    }
                }
             #   cat("1","\n")
                df<-data.frame("NX_code"= rep(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"])), length(output_total)), "GO_code"=paste(output_total) )

                xml_data_df<-rbind(xml_data_df,df)

            }else if(length(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"]))) ==1 & (length(unique(paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"]))) ==0)){
                
                
                subset_GO<-"NA"

                df<-data.frame("NX_code"= rep(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"])), length(subset_GO)), "GO_code"=paste(subset_GO) )
              #  cat("2","\n")
                xml_data_df<-rbind(xml_data_df,df)

            
            }

        }
        #### ultimo caso
        tmp<-xml_data_unlist[index_fin:length(xml_data_unlist)]
        tmp_df <- data.frame("name" = names(tmp), "value" = tmp)

            if(length(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"]))) ==1 & (length(unique(paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"]))) >=1)){

                subset_GO_tmp<-paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"])
                subset_KW<-subset_GO_tmp[startsWith(subset_GO_tmp, "KW")]
                if(length(subset_KW)!=0){
                    output_KW<-c()
                   for(KW in subset_KW){
                       #cat("KW ", KW,"\n")
                       sys.sleep(2)
                       file<-readLines(paste0("https://www.uniprot.org/keywords/",KW), warn=F)
                       GO_out<-str_extract(grep("QuickGO", file, value=T),'GO:[0-9]+')
                        if(length(GO_out)==1){
                            output_KW<-c(output_KW, GO_out)
                        }else if(length(GO_out) > 1){
                            #cat("KW distinto de 1",KW, "\n")
                        }else{
                            next
                        } 
                     
                   }

                }
                subset_GO<-subset_GO_tmp[startsWith(subset_GO_tmp, "GO:")]
                subset_GO<-unique(subset_GO)
                if(exists("output_KW")){
                        if(length(subset_GO)==0 & length(output_KW)==0){
                            output_total<-"NA"
                        }else if(length(subset_GO) >0 & length(output_KW) >0){
                            output_total<-c(subset_GO, output_KW)
                            output_total<-unique(output_total)
                        }else if(length(subset_GO) >0 & length(output_KW)==0){
                            output_total<-subset_GO
                        }else if(length(subset_GO)==0 & length(output_KW) >0){
                            output_total<-unique(output_KW)
                        } 
                }else{
                    if(length(subset_GO)==0){
                    
                        output_total<-"NA"
                    }else{
                        output_total<-subset_GO
                    }
                }
             #   cat("1","\n")
                df<-data.frame("NX_code"= rep(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"])), length(output_total)), "GO_code"=paste(output_total) )

                xml_data_df<-rbind(xml_data_df,df)

               # cat("3","\n")
            }else if(length(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"]))) ==1 & (length(unique(paste(tmp_df[tmp_df$name == "entry-list.entry.annotation-list.annotation-category.annotation.cv-term..attrs.accession","value"]))) ==0)){
                
                subset_GO<-"NA"
                #cat("4","\n")
                df<-data.frame("NX_code"= rep(unique(paste(tmp_df[tmp_df$name =="entry-list.entry..attrs.accession","value"])), length(subset_GO)), "GO_code"=paste(subset_GO) )

                xml_data_df<-rbind(xml_data_df,df)

            }
        

    return(xml_data_df)

}


data_f<-xmlParse(paste0(directory,xmlFile_in))
xml_data_raw <- xmlToList(data_f)
xml_data_unlist <- unlist(xml_data_raw)


Nextprot_output<-parsedataXML(xml_data_unlist)


name_out<-lapply(strsplit(paste(xmlFile_in),"\\."),"[",1)

write.table(Nextprot_output, file=paste0(directory,name_out,".tsv"), col.names=T, row.names=F, quote=F, sep="\t")

cat(system("date"),"\n")

kk<-paste0("DONE ", paste(xmlFile_in), "\n")
cat(kk,"\n")



