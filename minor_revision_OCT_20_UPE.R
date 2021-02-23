
setwd("/home/nostromo/data/03_Analysis/eguruce/20_UPs_Abr19")

load("EvaluationOfPR_sep20.RData")

### PR_CAFA_BP distintas de 1



noexact_match_CAFA_CCLEBP<-PE10_PRvsCAFA_CCLEBP_max[which(PE10_PRvsCAFA_CCLEBP_max$Dist.max <1),] #700/1304
noexact_match_CAFA_CCLECC<-PE10_PRvsCAFA_CCLECC_max[which(PE10_PRvsCAFA_CCLECC_max$Dist.max <1),] #680/1407
noexact_match_CAFA_CCLEMF<-PE10_PRvsCAFA_CCLEMF_max[which(PE10_PRvsCAFA_CCLEMF_max$Dist.max <1),] #781/1064


noexact_match_CAFA_GTEXBP<-PE10_PRvsCAFA_GTEXBP_max[which(PE10_PRvsCAFA_GTEXBP_max$Dist.max <1),] #648/1307
noexact_match_CAFA_GTEXCC<-PE10_PRvsCAFA_GTEXCC_max[which(PE10_PRvsCAFA_GTEXCC_max$Dist.max <1),] #677/1427
noexact_match_CAFA_GTEXMF<-PE10_PRvsCAFA_GTEXMF_max[which(PE10_PRvsCAFA_GTEXMF_max$Dist.max <1),] #717/1077

noexact_match_CAFA_TCGABP<-PE10_PRvsCAFA_TCGABP_max[which(PE10_PRvsCAFA_TCGABP_max$Dist.max <1),] #653/1075
noexact_match_CAFA_TCGACC<-PE10_PRvsCAFA_TCGACC_max[which(PE10_PRvsCAFA_TCGACC_max$Dist.max <1),] #588/1421
noexact_match_CAFA_TCGAMF<-PE10_PRvsCAFA_TCGAMF_max[which(PE10_PRvsCAFA_TCGAMF_max$Dist.max <1),] #653/1075


ALL_NOEXACT<-rbind(noexact_match_CAFA_CCLEBP,noexact_match_CAFA_CCLECC,noexact_match_CAFA_CCLEMF,noexact_match_CAFA_GTEXBP,noexact_match_CAFA_GTEXCC,noexact_match_CAFA_GTEXMF,noexact_match_CAFA_TCGABP,noexact_match_CAFA_TCGACC,noexact_match_CAFA_TCGAMF )

ALL_NOEXACT$Protein<-lapply(strsplit(paste(ALL_NOEXACT$Prot),"_GO"),"[",1)
length(unique(paste(ALL_NOEXACT$Protein))) # 1427




all_PE10_PR_CAFA<-rbind(PE10_PRvsCAFA_CCLEBP_max,PE10_PRvsCAFA_CCLECC_max, PE10_PRvsCAFA_CCLEMF_max,PE10_PRvsCAFA_GTEXBP_max, PE10_PRvsCAFA_GTEXCC_max,PE10_PRvsCAFA_GTEXMF_max, PE10_PRvsCAFA_TCGABP_max, PE10_PRvsCAFA_TCGACC_max, PE10_PRvsCAFA_TCGAMF_max)
all_PE10_PR_CAFA$Protein<-lapply(strsplit(paste(all_PE10_PR_CAFA$Prot),"_GO"),"[",1)
length(unique(paste(all_PE10_PR_CAFA$Protein))) #1560



#ir a nextprot y rescatar predicciones para esas proteinas



#### predicciones en Nextprot para estas proteinas. Sacaré también las MF que hay en Nextprot

all_PR_NX<-rbind(PE10_PRvsNx_CCLEBP_max,PE10_PRvsNx_CCLECC_max,PE10_PRvsNx_CCLEMF_max, PE10_PRvsNx_GTEXBP_max, PE10_PRvsNx_GTEXCC_max,PE10_PRvsNx_GTEXMF_max, PE10_PRvsNx_TCGABP_max, PE10_PRvsNx_TCGACC_max, PE10_PRvsNx_TCGAMF_max)

# ver las proteinas que no damos en el clavo en NX

all_PR_NX$Protein<-lapply(strsplit(paste(all_PR_NX$Prot),"_GO"),"[",1)

all_PR_NXselected<-all_PR_NX[which(all_PR_NX$Protein %in% ALL_NOEXACT$Protein),]
all_PR_NXselected_dist1<-all_PR_NXselected[which(all_PR_NXselected$Dist.max==1),]
length(unique(paste(all_PR_NXselected_dist1$Protein))) #1376



########################


VP & FP & FN & TPR & PPV & FNR & FDR & F1 & MÉTODO & REFERENCIA
867 & 48 & 118 & 0.88 & 0.94 & 0.12 & 0.05 & 0.91 & MACS & EXP1
804 & 483 & 181 & 0.81 & 0.62 & 0.18 & 0.37 & 0.70 & BAYES & EXP1
802 & 48 & 183 & 0.81 & 0.94 & 0.18 & 0.05 & 0.87 & JAMM & EXP1
953 & 332 & 32 & 0.96 & 0.74 & 0.03 & 0.25 & 0.83 & SICER & EXP1
881 & 199 & 104 & 0.89 & 0.81 & 0.10 & 0.18 & 0.85 & RANGER & EXP1
710 & 80 & 275 & 0.72 & 0.89 & 0.27 & 0.10 & 0.81 & ZCL & EXP1
1509 & 108 & 64 & 0.95 & 0.93 & 0.04 & 0.06 & 0.94 & MACS & EXP2
1476 & 513 & 97 & 0.93 & 0.74 & 0.06 & 0.25 & 0.82 & BAYES & EXP2
1420 & 114 & 153 & 0.90 & 0.92 & 0.09 & 0.07 & 0.91 & JAMM & EXP2
1562 & 394 & 11 & 0.99 & 0.79 & 0.01 & 0.20 & 0.88 & SICER & EXP2
1520 & 270 & 53 & 0.96 & 0.84 & 0.03 & 0.15 & 0.90 & RANGER & EXP2
1298 & 111 & 275 & 0.82 & 0.92 & 0.17 & 0.07 & 0.87 & ZCL & EXP2
1061 & 122 & 137 & 0.88 & 0.89 & 0.11 & 0.10 & 0.89 & MACS & EXP3
992 & 551 & 206 & 0.82 & 0.64 & 0.17 & 0.35 & 0.72 & BAYES & EXP3
943 & 124 & 255 & 0.78 & 0.88 & 0.21 & 0.11 & 0.83 & JAMM & EXP3
1158 & 388 & 40 & 0.96 & 0.74 & 0.03 & 0.25 & 0.84 & SICER & EXP3
1053 & 279 & 145 & 0.87 & 0.79 & 0.12 & 0.20 & 0.83 & RANGER & EXP3
885 & 118 & 313 & 0.73 & 0.88 & 0.26 & 0.11 & 0.80 & ZCL & EXP3