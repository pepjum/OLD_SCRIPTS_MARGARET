




data<-read.table("/home/margaret/data/pepe/GOannotation_Ago20.txt", header=T, sep="\t")

dir.create("/home/margaret/data/pepe/GOATOOLS/")

GOIDS<-paste(data$GOID)    # 47344
setwd("~/data/pepe/GOATOOLS_v2013/")
for(GO in GOIDS){
	y<-paste0("python3 ~/.local/bin/wr_hier.py --dag /home/margaret/data/pepe/gene_ontology_2013.obo -f ", GO)
	system(y)
}

### load files
library(stringr)

files<-list.files(pattern=".txt")


df_output<-data.frame()

for (file in files){
	cat(file,"\n")
	file_out<-read.table(file, sep="\t")

	GOname<-basename(file)
	number<-str_extract(GOname,"[0-9]+")
	GO_number<-paste0("GO:",number)


	selectedRows <- file_out[grep(GO_number, file_out$V1), ]
	if(length(selectedRows)==1){
		level<-str_extract(selectedRows,"L[0-9]+")
		output<-data.frame("GO"=GO_number,"level"=level)
		df_output<-rbind(df_output, output)
	}else{
		levels<-str_extract(selectedRows,"L[0-9]+")
		levels_sorted<-sort(levels, decreasing=T)
		level<-levels_sorted[1]
		output<-data.frame("GO"=GO_number,"level"=level)
		df_output<-rbind(df_output, output)
	}
}


write.table(df_output, file="~/data/pepe/GO_LEVELS_2013.txt", col.names=F, row.names=F, sep="\t", quote=F)


#### collect datas

GO_number<-str_extract(GOIDS,"[0-9]+")
namefile<-paste0("hier_GO",GO_number,".txt")

not_exist<-c()
files<-paste0("/home/margaret/data/pepe/GOATOOLS/",namefile)
#files_old_DB<-paste0("/home/margaret/data/pepe/GOATOOLS_v2013/", namefile)

df_output<-data.frame()
for(file in files){
	
	if(file.exists(file)){
		info<-file.info(file)
		if(info$size >0){
			cat(file,"\n")
			file_out<-readLines(file)
			file_out<-as.data.frame(file_out)
			GOname<-basename(file)
			number<-str_extract(GOname,"[0-9]+")
			GO_number<-paste0("GO:",number)

			colnames(file_out)<-"V1"
			selectedRows <- file_out[grep(GO_number, file_out$V1), ]
				if(length(selectedRows)==1){
					level<-str_extract(selectedRows,"L[0-9]+")
					output<-data.frame("GO"=GO_number,"level"=level)
					df_output<-rbind(df_output, output)
				}else{
					levels<-str_extract(selectedRows,"L[0-9]+")
					levels_sorted<-sort(levels, decreasing=T)
					level<-levels_sorted[1]
					output<-data.frame("GO"=GO_number,"level"=level)
					df_output<-rbind(df_output, output)
				}
		
		}else{
			

			nf<-basename(file)
			queryfile<-paste0("/home/margaret/data/pepe/GOATOOLS_v2013/", nf)
			info<-file.info(queryfile)
			if(info$size >0){
				cat(queryfile,"\n")

				if(file.exists(queryfile)){
					file_out<-readLines(queryfile)
					file_out<-as.data.frame(file_out)
					GOname<-basename(file)
					number<-str_extract(GOname,"[0-9]+")
					GO_number<-paste0("GO:",number)

					colnames(file_out)<-"V1"

					selectedRows <- file_out[grep(GO_number, file_out$V1), ]
						if(length(selectedRows)==1){
							level<-str_extract(selectedRows,"L[0-9]+")
							output<-data.frame("GO"=GO_number,"level"=level)
							df_output<-rbind(df_output, output)
						}else{
							levels<-str_extract(selectedRows,"L[0-9]+")
							levels_sorted<-sort(levels, decreasing=T)
							level<-levels_sorted[1]
							output<-data.frame("GO"=GO_number,"level"=level)
							df_output<-rbind(df_output, output)
						}

				}else{
					not_exist<-c(not_exist,queryfile)
				}
			}else{
				not_exist<-c(not_exist, queryfile)
			}
		}
	}else{
		
		nf<-basename(files)
		queryfile<-paste0("/home/margaret/data/pepe/GOATOOLS_v2013/", nf)
		info<-file.info(queryfile)
		if(info$size >0){
			cat(queryfile,"\n")

			if(file.exists(queryfile)){
				file_out<-readLines(queryfile)
				file_out<-as.data.frame(file_out)
				GOname<-basename(file)
				number<-str_extract(GOname,"[0-9]+")
				GO_number<-paste0("GO:",number)

				colnames(file_out)<-"V1"

				selectedRows <- file_out[grep(GO_number, file_out$V1), ]
					if(length(selectedRows)==1){
						level<-str_extract(selectedRows,"L[0-9]+")
						output<-data.frame("GO"=GO_number,"level"=level)
						df_output<-rbind(df_output, output)
					}else{
						levels<-str_extract(selectedRows,"L[0-9]+")
						levels_sorted<-sort(levels, decreasing=T)
						level<-levels_sorted[1]
						output<-data.frame("GO"=GO_number,"level"=level)
						df_output<-rbind(df_output, output)
					}

			}else{
				not_exist<-c(not_exist,queryfile)
			}
		}else{
			not_exist<-c(not_exist, queryfile)
		}

	}
}


faltan<-setdiff(GOIDS, df_output$GO)

recuperados<-data.frame()
counter<-0
for(GO in faltan){
	counter<-counter+1
	cat(counter,"/",length(faltan),"\n")
    file<-readLines(paste0('http://golr-aux.geneontology.io/solr/select?fq=document_category:"ontology_class"&q=*:*&fq=id:','"',GO,'"','&wt=json'), warn=F)
	filepart<-lapply(strsplit(paste(file),"replaced_by"),"[",2)
	if(!is.na(filepart)){
    	GO_out<-str_extract(filepart,'GO:[0-9]+')
		reemplazados<-data.frame("old"=GO,"new"=GO_out)
		recuperados<-rbind(recuperados,reemplazados)
	}else{

		filepart<-lapply(strsplit(paste(file),"consider"),"[",2)
		if(!is.na(filepart)){
			GO_out<-str_extract(filepart,'GO:[0-9]+')[1]

			no_reemplazados<-data.frame("old"=GO,"new"=GO_out)
			recuperados<-rbind(recuperados,no_reemplazados)
		}else{
			no_reemplazados<-data.frame("old"=GO, "new"="OBSOLETE")
			recuperados<-rbind(recuperados, no_reemplazados)
		}
	}
}

recuperados[which(recuperados$old=="GO:0100022"),]$new<-"OBSOLETE"
obsoletos<-recuperados[which(recuperados$new=="OBSOLETE"),]
colnames(obsoletos)<-c("GO","level")
df_output<-rbind(df_output,obsoletos)

no_obsoletos<-recuperados[which(recuperados$new!="OBSOLETE"),]

no_obsoletos_level<-merge(no_obsoletos,df_output, by.x="new", by.y="GO")

no_obsoletos_level_sel<-no_obsoletos_level[,c("old","level")]
names(no_obsoletos_level_sel)<-c("GO","level")
df_output<-rbind(df_output,no_obsoletos_level_sel)

setdiff(GOIDS,df_output$GO)   #"GO:0008557"
uno<-data.frame(old="GO:0008557", new="GO:0140326")

df_output<-rbind(df_output,c("GO:0008557","L04"))

