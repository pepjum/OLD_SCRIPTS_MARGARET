


fileraw<-"/home/margaret/data/pepe/pruebas_simulados/MSMSPredict_2.txt"

con <- file(fileraw, open="r")
simulated_peptides <- readLines(con)
close(con)

options(scipen=999)

# grepl('[A-Z]+',simulated_peptides[i] )
# grepl('[0-9]+',simulated_peptides[i] )
#
# tmp<-c()
tmp<-simulated_peptides[-c(1)]
tmp<-lapply(tmp, FUN=function(x) x=paste0(x,"####")) #meto un fin de linea reconocible
for(i in 2:length(tmp)){
    if(grepl("[A-Z]+", tmp[i]) & grepl("[0-9]+", tmp[i-1])){
        tmp[i-1]<-paste0(tmp[i-1],"@")
    }
}

simulated_list<-paste(unlist(tmp), collapse="")
tmp_l<-unlist(strsplit(paste(simulated_list),"@"))
tmp_l <- unlist(lapply(tmp_l, FUN = function(x) x[x != ""]))
tmp_l2<-unlist(strsplit(paste(tmp_l),"####"))

peptides<-c()
charge<-c()
description<-c()
method<-c()
Collision_energy<-c()
Reaction_time<-c()
isolation_width<-c()
instrument_model<-c()
resolution_of_peptide<-c()
scan_range<-c()
total_number_of_mass_intensity<-c()
Mass_intensity_pairs<-c()

for(i in 1:length(tmp_l)){
    string<-unlist(strsplit(tmp_l[i],"####"))
    peptides<-c(peptides,string[1])
    charge<-c(charge,string[3])
    method<-c(method,string[4])
    Collision_energy<-c(Collision_energy,string[5])
    Reaction_time<-c(Reaction_time,string[6])
    isolation_width<-c(isolation_width, string[7])
    instrument_model<-c(instrument_model,string[8])
    resolution_of_peptide<-c(resolution_of_peptide, string[9])
    scan_range<-c(scan_range,string[10])
    total_number_of_mass_intensity<-c(total_number_of_mass_intensity,string[11])
    Mass_intensity_pairs<-c(Mass_intensity_pairs,paste(string[12:length(string)],collapse="##"))
}

mass_df<-data.frame("Peptide"=unique(peptides))
# mass_df$mass<-as.numeric(c("2776.2658","4592.1483","3443.7667","1078.4952","1030.5680","1283.63","840.4106","1853.9214","696.3749","1214.5946","778.3413","1608.7798","687.3784","909.4312","3877.6205","2068.1760","751.3733","705.3171","890.4214","1385.6808","1340.6804","1579.8478","1058.5728", "589.2940","602.2892","779.4158", "635.3107", "720.2682"))
mass_df$mass<-c(570.3059,1011.57996,821.42255)
save(mass_df, file="/home/margaret/data/pepe/mass_correlation_peptides_TP53_simulados_unicos.rda")

#charge2<-charge[c(1,10,15,17,19,21,22,25,26,29,30,33,34,35,40,44,45,46,47,50,52,55,57,58,59,60,61)]
#peptides2<-peptides[c(1,10,15,17,19,21,22,25,26,29,30,33,34,35,40,44,45,46,47,50,52,55,57,58,59,60,61)]


tmp_output<-data.frame("peptides"=paste(peptides),"charge"=paste(charge))

tmp_2_output<-merge(tmp_output, mass_df, by.x="peptides", by.y="Peptide", all.x=T)
#tmp_2_output$mass_cor<-as.numeric(paste(tmp_2_output$mass))/as.numeric(paste(tmp_2_output$charge))+1
tmp_2_output$intensities<-paste(Mass_intensity_pairs)
tmp_2_output$scan<-seq(1:nrow(tmp_2_output))

list_dataframes<-list()
monta_df<-function(objeto){
    a<-objeto
    b<-unlist(strsplit(a,"##"))
    m<-sapply(strsplit(b,"\t"),"[",1)
    z<-sapply(strsplit(b,"\t"),"[",2)
    c<-data.frame("m"=m,"z"=z)
}

#Mass_intensity_pairs2<-Mass_intensity_pairs[c(1,10,15,17,19,21,22,25,26,29,30,33,34,35,40,44,45,46,47,50,52,55,57,58,59,60,61)]


for (i in 1:length(Mass_intensity_pairs)){
    objeto<-Mass_intensity_pairs[i]
    df_out<-monta_df(objeto)
    list_dataframes[[i]]<-df_out
}
output_list<-c()
for(i in 1:nrow(tmp_2_output)){
    spec<-tmp_2_output[i,]
    a<-"BEGIN IONS"
    output_list<-c(output_list,a)
    Title<-paste0("TITLE=",spec$scan)
    output_list<-c(output_list,Title)
    Pepmass<-paste0("PEPMASS=", spec$mass/2)
    output_list<-c(output_list,Pepmass)
    Charge<-paste0("CHARGE=",spec$charge)
    output_list<-c(output_list,Charge)
    Seqs<-paste0("SEQ=",spec$peptides)
    output_list<-c(output_list,Seqs)
    peaks<-as.data.frame(list_dataframes[i])
    for (k in 1:nrow(peaks)){
        lines<-peaks[k,]
        line<-paste(lines$m,lines$z, sep="\t")
        output_list<-c(output_list,line)
    }
    b<-"END IONS\n"
    output_list<-c(output_list,b)
}

writeLines(unlist(lapply(output_list, paste, collapse="\n")), "/home/margaret/data/pepe/TP53_simulado_con_masa_massanalyzer.mgf")
