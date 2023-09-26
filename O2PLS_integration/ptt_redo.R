library(readxl)
library(dplyr)
library(writexl)

Protein_All_Results <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/Protein_All_Results.xlsx")
prt_FDR_adj = Protein_All_Results[Protein_All_Results$adj.P.Val < 0.05,]
write_xlsx(prt_FDR_adj,"prt_FDR_adj.xlsx")

####without using o2pls ###

###rnaseq filtering
##to separate proein directions
#get the direction of joint loading values of the final 
ptt_jnt1_FDR_adj <- read_excel("~/BINF/multi-omics/o2pls_all_data/post transcriptional regulation/ptt_jnt1_FDR_adj.xlsx")
ptt_jnt2_FDR_adj <- read_excel("~/BINF/multi-omics/o2pls_all_data/post transcriptional regulation/ptt_jnt2_FDR_adj.xlsx")

#load FC of all prots

Protein_All_Results <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/Protein_All_Results.xlsx", 
                                     sheet = "avg_FC")

ptt_jnt1_FDR_adj_FC = inner_join(ptt_jnt1_FDR_adj,Protein_All_Results, by = "ID")
write_xlsx(ptt_jnt1_FDR_adj_FC, "ptt_jnt1_FDR_adj_FC.xlsx")

ptt_jnt2_FDR_adj_FC = inner_join(ptt_jnt2_FDR_adj,Protein_All_Results, by = "ID")
write_xlsx(ptt_jnt2_FDR_adj_FC, "ptt_jnt2_FDR_adj_FC.xlsx")

#separate positive and negative FCs

ptt_jnt1_FDR_adj_FC_neg <- ptt_jnt1_FDR_adj_FC[ptt_jnt1_FDR_adj_FC$logFC < 0,]
write_xlsx(ptt_jnt1_FDR_adj_FC_neg, "ptt_jnt1_FDR_adj_FC_neg.xlsx")
ptt_jnt2_FDR_adj_FC_neg <- ptt_jnt2_FDR_adj_FC[ptt_jnt2_FDR_adj_FC$logFC < 0,]
write_xlsx(ptt_jnt2_FDR_adj_FC_neg,"ptt_jnt2_FDR_adj_FC_neg.xlsx")


ptt_jnt1_FDR_adj_FC_pos <- ptt_jnt1_FDR_adj_FC[ptt_jnt1_FDR_adj_FC$logFC > 0,]
write_xlsx(ptt_jnt1_FDR_adj_FC_pos, "ptt_jnt1_FDR_adj_FC_pos.xlsx")
ptt_jnt2_FDR_adj_FC_pos <- ptt_jnt2_FDR_adj_FC[ptt_jnt2_FDR_adj_FC$logFC > 0,]
write_xlsx(ptt_jnt2_FDR_adj_FC_pos,"ptt_jnt2_FDR_adj_FC_pos.xlsx")

##filter with surfaceome
surfaceome <- read_excel("~/BINF/secretomics/references/table_S3_surfaceome.xlsx", 
                                   sheet = "in silico surfaceome only")
#change column name in surfaceome

names(surfaceome)[names(surfaceome) == "UniProt accession"] <- "ID"

ptt_jnt1_FDR_adj_FC_pos_surf = inner_join(ptt_jnt1_FDR_adj_FC_pos, surfaceome, by = "ID")
write_xlsx(ptt_jnt1_FDR_adj_FC_pos_surf, "ptt_jnt1_FDR_adj_FC_pos_surf.xlsx")
ptt_jnt2_FDR_adj_FC_pos_surf = inner_join(ptt_jnt2_FDR_adj_FC_pos, surfaceome, by = "ID")
write_xlsx(ptt_jnt2_FDR_adj_FC_pos_surf, "ptt_jnt2_FDR_adj_FC_pos_surf.xlsx")

common_substrates_jnt12 = inner_join(ptt_jnt1_FDR_adj_FC_pos_surf, ptt_jnt2_FDR_adj_FC_pos_surf, by = "ID")
common_substrates_jnt12_names = common_substrates_jnt12[,c(1,2)]
write_xlsx(common_substrates_jnt12_names, "common_substrates_jnt12_names.xlsx")

#how many of these 93 can be found in the original 179 potential substrates?
original = read.table(file ="clipboard", header = TRUE)
from_original = ptt_jnt2_FDR_adj_FC_pos_surf[ptt_jnt2_FDR_adj_FC_pos_surf$Symbol.x %in% original$gene,]
nn= anti_join(ptt_jnt2_FDR_adj_FC_pos_surf, from_original)

#find how many common exist

pos_no_common = inner_join(ptt_jnt1_FDR_adj_FC_pos, ptt_jnt2_FDR_adj_FC_pos, by = "ID") #1203
pos_union = full_join(ptt_jnt1_FDR_adj_FC_pos, ptt_jnt2_FDR_adj_FC_pos)
pos_union_names = pos_union[,c(1,2,4)]
write_xlsx(pos_union_names,"pos_union_names.xlsx")

#negative FC
ptt_jnt1_FDR_adj_FC_neg <- read_excel("ptt_jnt1_FDR_adj_FC_neg.xlsx")
ptt_jnt2_FDR_adj_FC_neg <- read_excel("ptt_jnt2_FDR_adj_FC_neg.xlsx")
#neg_commons = anti_join(ptt_jnt1_FDR_adj_FC_neg,ptt_jnt2_FDR_adj_FC_neg, by="ID" )
neg_union = full_join(ptt_jnt1_FDR_adj_FC_neg, ptt_jnt2_FDR_adj_FC_neg, by="ID")
neg_union_names = neg_union2[,c(1,2,4)]
write_xlsx(neg_union_names, "neg_union_names.xlsx")
neg_union2 = full_join(ptt_jnt1_FDR_adj_FC_neg, ptt_jnt2_FDR_adj_FC_neg)

#different fold change cut-offs for PTTs

#neg_PTTs
neg_un_1.2 = neg_union_names[neg_union_names$logFC < -1.2,]
write_xlsx(neg_un_1.2,"neg_un_1.2.xlsx")

neg_un_1 = neg_union_names[neg_union_names$logFC < -1,]
write_xlsx(neg_un_1,"neg_un_1.xlsx")

neg_un_1.5 = neg_union_names[neg_union_names$logFC < -1.5,]
write_xlsx(neg_un_1.5,"neg_un_1.5.xlsx")

neg_un_0.5 = neg_union_names[neg_union_names$logFC < -0.5,]
write_xlsx(neg_un_0.5,"neg_un_0.5.xlsx")

#pos_PTT
pos_un_1.2 = pos_union_names[pos_union_names$logFC > 1.2,]
write_xlsx(pos_un_1.2,"pos_un_1.2.xlsx")

pos_un_1.5 = pos_union_names[pos_union_names$logFC > 1.5,]
write_xlsx(pos_un_1.5,"pos_un_1.5.xlsx")


pos_un_1 = pos_union_names[pos_union_names$logFC > 1,]
write_xlsx(pos_un_1,"pos_un_1.xlsx")


pos_un_0.5 = pos_union_names[pos_union_names$logFC > 0.5,]
write_xlsx(pos_un_0.5,"pos_un_0.5.xlsx")

