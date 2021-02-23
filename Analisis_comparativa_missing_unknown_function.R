#### Cruzar listas de nextprot con funcion desconocida con las missing detectadas en embrionarias



missing<-read.table("~/data/pepe/42_EMBRIO_FUNCIONS_JUL20/all_missing_proteins_detected.txt")
nextprot_unknown_function<-read.table("~/data/pepe/42_EMBRIO_FUNCIONS_JUL20/nextprot-unknown_function.txt")


missing_2_or_more_peptides<-missing[which(missing$V2 >1),]  # 13 missing

missing_unknown_function_2_or_more_peptides<-missing_2_or_more_peptides[which(missing_2_or_more_peptides$V1 %in% paste(nextprot_unknown_function$V1)),]   # 3 missing

write.table(missing_unknown_function_2_or_more_peptides, file="~/data/pepe/42_EMBRIO_FUNCIONS_JUL20/missing_2_or_more_peptides_unknown_function.txt", col.names=F, row.names=F, quote=F, sep="\t")

one_hit<-missing[which(missing$V2 ==1),]

one_hit_unknown_function_1_peptide<-one_hit[which(one_hit$V1 %in% paste(nextprot_unknown_function$V1)),]   # 27 one hit with unknown function

write.table(one_hit_unknown_function_1_peptide, file="~/data/pepe/42_EMBRIO_FUNCIONS_JUL20/one_hit_wonders_unknown_function.txt", col.names=F, row.names=F, quote=F, sep="\t")
