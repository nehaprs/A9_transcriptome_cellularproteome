library(dplyr)
library(writexl)
library(readxl)
#load lvs of joint component 1

rnaseq_cp_all_joint <- read_csv("~/BINF/multi-omics/o2pls_all_data/rnaseq_cp/results/common_rnaseq_cp_all_joint.csv")
row =  row = rnaseq_cp_all_joint[1,]
g = row$V1gene
row[,4]
p = row$V1cp
if(sign(g) == sign(p)){print(1)} 
  
##joint component 1
#find those proteins who have lv 1 in proteomics have the opposite sign of that in transcriptomics

#create new dataset
rnaseq_cp_joint1 = rnaseq_cp_all_joint[,c(2,3,4,6)]
jnt1_opposite_signs = data.frame(gene = character(), ID = character(), jnt1gene = integer(), jnt1cp = integer())

for(i in 1:6262){
  row = rnaseq_cp_joint1[i,]
  g = row$V1gene
  p = row$V1cp
  if(sign(g) != sign(p)){
    jnt1_opposite_signs = rbind(jnt1_opposite_signs, row) 
  } 
}

write_xlsx(jnt1_opposite_signs, "jnt1_opposite_signs")

##do the same for joint component 2
#find those proteins who have lv 2 in proteomics have the opposite sign of that in transcriptomics

#create new dataset
rnaseq_cp_joint2 = rnaseq_cp_all_joint[,c(2,3,5,7)]
jnt2_opposite_signs = data.frame(gene = character(), ID = character(), jnt1gene = integer(), jnt1cp = integer())

for(i in 1:6262){
  row = rnaseq_cp_joint2[i,]
  g = row$V2gene
  p = row$V2cp
  if(sign(g) != sign(p)){
    jnt2_opposite_signs = rbind(jnt2_opposite_signs, row) 
  } 
}

write_xlsx(jnt2_opposite_signs, "jnt2_opposite_signs")

##rescuing more proteins using Pearson

##find the prots that are in all but not opposite sign dataset

to_rescue_jnt1 = anti_join(rnaseq_cp_joint1,jnt1_opposite_signs, by = "ID")
to_rescue_jnt2 = anti_join(rnaseq_cp_joint2,jnt2_opposite_signs, by = "ID")
write_xlsx(to_rescue_jnt1,"to_rescue_jnt1.xlsx")
write_xlsx(to_rescue_jnt2,"to_rescue_jnt2.xlsx")

#find subset of original data with these as rows

##load original data with FC

mRNA_All_Results <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx", 
                                  sheet = "fc")

Protein_All_Results <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/Protein_All_Results.xlsx", 
                                     sheet = "fc")

mrna_fc = mRNA_All_Results[,c(1,14,15,16,17,18,19)]
prot_fc = Protein_All_Results[,c(1,2,13,14,15,16,17)]
mrna_fc = mrna_fc[,-7]
#select those rows to rescue
#for jc1

mrna_fc_toRescue1 = mrna_fc[mrna_fc$gene %in% to_rescue_jnt1$gene,]
prot_fc_toRescue1 = prot_fc[prot_fc$ID %in% to_rescue_jnt1$ID,]

#find pearson for these protein combos
#select first 5 rows from mrna 


#standardize the fold change values
fcrna_std = mrna_fc_toRescue1[,-c(1)] %>% mutate_all(~(scale(.) %>% as.vector))
fcprot_std = prot_fc_toRescue1[,-c(1:2)] %>% mutate_all(~(scale(.) %>% as.vector))

#extract names
fcrna_names = mrna_fc_toRescue1[,1]
fcprot_names = prot_fc_toRescue1[,c(1,2)]

#combine names and standardized values
fcrna_std_named = cbind(fcrna_names, fcrna_std)
fcprot_std_named = cbind(fcprot_names, fcprot_std)

#combine the standardized values to form 1 dataframe with standardized values
library(dplyr)
fcdata_std = inner_join(fcrna_std_named, fcprot_std_named, by ="gene")

# declaring an empty data frame to populate gene name, pearson, p-value
ToRescue_jnt1_pearson = data.frame(gene = character(), PearsonCoeff = numeric(), p_value = numeric())

#remove missing values before doing pearson
fcdata_std <- fcdata_std[complete.cases(fcdata_std), ]
#find pearson's correlation gene by gene

for(i in 1: nrow(fcdata_std))
{rowi = fcdata_std[i,]
ToRescue_jnt1_pearson[i,1] = rowi[1]
fcrnai = unlist(as.vector(rowi[,2:6]))
fcproti = unlist( as.vector(rowi[,8:12]))
pearsi = cor.test(fcrnai, fcproti, method = "pearson")
ToRescue_jnt1_pearson[i,2] = pearsi$estimate
ToRescue_jnt1_pearson[i,3] = pearsi$p.value

} 
write_xlsx(ToRescue_jnt1_pearson, "ToRescue_jnt1_pearson")

#####################

#repeat for joint component 2
mrna_fc_toRescue2 = mrna_fc[mrna_fc$gene %in% to_rescue_jnt2$gene,]
prot_fc_toRescue2 = prot_fc[prot_fc$ID %in% to_rescue_jnt2$ID,]

#find pearson for these protein combos
#select first 5 rows from mrna 


#standardize the fold change values
fcrna_std2 = mrna_fc_toRescue2[,-c(1)] %>% mutate_all(~(scale(.) %>% as.vector))
fcprot_std2 = prot_fc_toRescue2[,-c(1:2)] %>% mutate_all(~(scale(.) %>% as.vector))

#extract names
fcrna_names2 = mrna_fc_toRescue2[,1]
fcprot_names2 = prot_fc_toRescue2[,c(1,2)]

#combine names and standardized values
fcrna_std_named2 = cbind(fcrna_names2, fcrna_std2)
fcprot_std_named2 = cbind(fcprot_names2, fcprot_std2)

#combine the standardized values to form 1 dataframe with standardized values
library(dplyr)
fcdata_std2 = inner_join(fcrna_std_named2, fcprot_std_named2, by ="gene")

# declaring an empty data frame to populate gene name, pearson, p-value
ToRescue_jnt2_pearson = data.frame(gene = character(), PearsonCoeff = numeric(), p_value = numeric())

#remove missing values before doing pearson
fcdata_std2 <- fcdata_std2[complete.cases(fcdata_std2), ]
#find pearson's correlation gene by gene

for(i in 1: nrow(fcdata_std2))
{rowi = fcdata_std2[i,]
ToRescue_jnt2_pearson[i,1] = rowi[1]
fcrnai = unlist(as.vector(rowi[,2:6]))
fcproti = unlist( as.vector(rowi[,8:12]))
pearsi = cor.test(fcrnai, fcproti, method = "pearson")
ToRescue_jnt2_pearson[i,2] = pearsi$estimate
ToRescue_jnt2_pearson[i,3] = pearsi$p.value

} 
write_xlsx(ToRescue_jnt2_pearson, "ToRescue_jnt2_pearson.xlsx")

#########################################
#how much do the genes overlap between the components

combined_list_of_pttjnt1 <- read_excel("combined list of ptt.xlsx", 
                                          sheet = "jnt1")

combined_list_of_pttjnt2 <- read_excel("combined list of ptt.xlsx", 
                                       sheet = "jnt2")

combined_list_of_ptt_both = data.frame(intersect(combined_list_of_pttjnt1$gene ,combined_list_of_pttjnt2$gene))
write_xlsx(combined_list_of_ptt_both, "combined_list_of_ptt_both.xlsx")

##########################################

# of the combined list, how many have significant protein change opposite to adm9

prt_FDR_adj <- read_excel("redo/prt_FDR_adj.xlsx")
prtnames_FDR = prt_FDR_adj[,c(1,3)]

combined_list_of_ptt_jnt1 <- read_excel("combined list of ptt.xlsx", 
                                   sheet = "jnt1")
ptt_jnt1_FDR_adj = prtnames_FDR[prtnames_FDR$Symbol %in% combined_list_of_ptt_jnt1$gene,]

write_xlsx(ptt_jnt1_FDR_adj,"ptt_jnt1_FDR_adj.xlsx")
combined_list_of_ptt_jnt2 <- read_excel("combined list of ptt.xlsx", 
                                        sheet = "jnt2")

ptt_jnt2_FDR_adj = prtnames_FDR[prtnames_FDR$Symbol %in% combined_list_of_ptt_jnt2$gene,]
write_xlsx(ptt_jnt2_FDR_adj,"ptt_jnt2_FDR_adj.xlsx")
