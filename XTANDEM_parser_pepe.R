#!/usr/bin/Rscript

args=(commandArgs(TRUE))

xmlFile <- args[1]

cat(system("date"))
cat(xmlFile)
cat("\n")

require(XML)

data <- xmlParse(xmlFile)
xml_data_raw <- xmlToList(data)
save(xml_data_raw, file = paste(gsub(".xml","",xmlFile), "_raw.rda", sep = ""))
xml_data_unlist <- unlist(xml_data_raw)


parsedataXML<-function(xml_data_unlist){

	indexes <- which(names(xml_data_unlist) == "group..attrs.act")

	xml_data_df <- data.frame("Protein"=as.character(),
	"PeptideSeq"=as.character(),
	"scan"=as.character(),
	"PrecursorMZ"=as.character(),
	"charge"=as.character(),
	"expect"=as.character(),
	"index_mgf"=as.character(),
	"modification"=as.character(),
	"start"=as.character(),
	"end"=as.character(),
	"hyperscore"=as.character(),
	"missed_cleavages"=as.character())

	for (i in 1:length(indexes)) {

		if (i == 1) {
			index <- 1
			index_fin <- indexes[i] - 1
		}

		else {
			index <- indexes[i-1]
			index_fin <- indexes[i] - 1
		}

		tmp <- xml_data_unlist[index:index_fin]
		tmp_df <- data.frame("name" = names(tmp), "value" = tmp)

		# SI todo ES UNO rellenamos el dataframe del tiron
		if ((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) == 1) & ((length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"])))) == 1)))){


			 df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]),
				 "PeptideSeq"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])),
				 "scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
				 "PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
				 "charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
				 "expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
				 "index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
				 "modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"),
				 "start"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"])),
				 "end"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"])),
				 "hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"])), "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"]))

			 xml_data_df <- rbind(xml_data_df,df)

		}else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) > 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) == 1) &	(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) > 1) &
		(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"]))) > 1)){

			cat("loop 2","\n")
			cat(index, index_fin,"\n")
			break


			for (l in 1:length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]))){

				df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"][l]),
					"PeptideSeq"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])),
					"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
					"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
					"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
					"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
					"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
					"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"),
					"start"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][l]),
					"end"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][l]),
					"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"])), "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][l]))

				xml_data_df <- rbind(xml_data_df,df)

				}

		}
		# si tengo un protein id pero más de una secuencia
		else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) > 1) &	(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) > 1) &
		(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"]))) > 1)){

			cat("loop 3","\n")
			cat(index, index_fin,"\n")
			break


		 		df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]),
					"Peptideseq"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"][j]),
					"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
					"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
					"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
					"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
					"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
					"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"),
					"start"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][j]),
					"end"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][j]),
					"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"])), "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][l]))

		 		xml_data_df <- rbind(xml_data_df,df)
		 		}
		}
		#si tengo mas de un protein id y una secuencia
		else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) > 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end"
		,"value"]))) == 1)){

			cat("loop 4","\n")
			cat(index, index_fin,"\n")
			break


			for (j in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))){

				df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]),
					"PeptideSeq"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])),
					"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
					"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
					"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
					"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
					"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
					"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"), "start"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"])),
					"end"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"])),
					"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"])), "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][l]))

					xml_data_df <- rbind(xml_data_df,df)
				}
	    }
		else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) == 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))) > 1) & (length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end"
		,"value"]))) > 1)){

			cat("loop 5","\n")
			cat(index, index_fin,"\n")
			break


			for(m in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"]))){

				df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]),
					"PeptideSeq"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])),
					"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
					"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
					"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
					"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
					"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
					"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"), "start"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][m])),
					"end"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][m])),
					"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"][m])), "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][m]))

				xml_data_df <- rbind(xml_data_df,df)

			}
		}

		# si tengo más de un protein id y más de una secuencia
		else if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"])) > 1) &	(length(unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))) > 1)){
			# si esto van a uno a uno
			cat("loop 6","\n")
			cat(index, index_fin,"\n")
			break

			if((length(paste(tmp_df[tmp_df$name == "group.protein.note.text","value"]))) == (length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"])))){

					for (j in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))){
					#cat("\n")
					df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"][j]),
						"PeptideSeq"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"][j]),
						"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
						"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
						"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
						"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
						"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
						"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"), "start"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][j]),
						"end"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][j]),
						"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"][j])), "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][j]))


						xml_data_df <- rbind(xml_data_df,df)
				}
			}
			#si van mas de una secuencia por protein_id

			else{
				for (k in 1:length(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"]))){

					#cat("\n")
					df<- data.frame("Protein"=paste(tmp_df[tmp_df$name == "group.protein.note.text","value"][as.numeric(unlist(strsplit(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.id","value"][k]),"[.]"))[2])]),
						"PeptideSeq"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.seq","value"][k]),
						"scan"=paste(tmp_df[tmp_df$name == "group.group.note.text","value"]),
						"PrecursorMZ"=paste(tmp_df[tmp_df$name == "group..attrs.mh","value"]),
						"charge"=paste(tmp_df[tmp_df$name == "group..attrs.z","value"]),
						"expect"=paste(tmp_df[tmp_df$name == "group..attrs.expect","value"]),
						"index_mgf"=paste(tmp_df[tmp_df$name == "group..attrs.id","value"]),
						"modification"=paste(paste(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.type","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.at","value"]),paste(tmp_df[tmp_df$name == "group.protein.peptide.domain.aa.modified","value"]),sep=","),collapse=";"), "start"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.start","value"][k]),
						"end"=paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.end","value"][k]),
						"hyperscore"=unique(paste(tmp_df[tmp_df$name == "group.protein.peptide.domain..attrs.hyperscore","value"][j])), "missed_cleavages"=paste(tmp_df[tmp_df$name=="group.protein.peptide.domain..attrs.missed_cleavages","value"][j]))



					xml_data_df <- rbind(xml_data_df,df)
				}
			}
		}
	}
}
kk<-parsedataXML(xml_data_unlist)
#save(xml_data_df, file = paste(gsub(".xml","",xmlFile), ".rda", sep = ""))

write.table(xml_data_df, file = paste(gsub(".xml","",xmlFile), ".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
cat(system("date"))



### mh -> group..attrs.mh
### charge -> group..attrs.z
# expect -> group..attrs.expect
# expectro -> group..attrs.id
