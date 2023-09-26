neha = read.table(file="clipboard", sep = "\t", header = T) 
congyu = read.table(file="clipboard", sep = "\t", header = T) 

library(dplyr)
intersection = inner_join(neha, congyu, by ="gene")
write.table(intersection,"congyu61overlap")
commonsJn1 = read.table(file="clipboard", sep = "\t", header = T, row.names = 1)
commonsJn2 = read.table(file="clipboard", sep = "\t", header = T, row.names = 1)
commonOfcommons = inner_join(commonsJn1, commonsJn2, by = 'genes')


surfaceome = read.table(file="clipboard", sep = "\t", header = T)
cpOrth =  read.table(file="clipboard", sep = "\t", header = T)
surface_in_cp_orth = inner_join(cpOrth, surfaceome, by ='ID')
write.table(surface_in_cp_orth,"surface_in_cp_orth")

high_cp_orth = read.table(file="clipboard", sep = "\t", header = T)
are_highcporthinner = inner_join(high_cp_orth, congyu, by = 'gene')

datajnt =  read.table(file="clipboard", sep = "\t", header = T, row.names = 1)
dataorth =  read.table(file="clipboard", sep = "\t", header = T, row.names = 1)
jntorth = inner_join(datajnt, dataorth)

common61 = read.table(file="clipboard", sep = "\t", header = T, row.names = 1)
jntcomp2vals = read.table(file="clipboard", sep = "\t", header = T, row.names = 1)
common61_all_vals = inner_join(common61, jntcomp2vals)
write.table(common61_all_vals, "Common61genes_withLoadingVals")

not_common61 = anti_join(congyu, common61)

# transcriptional regulation
#load all jnt1 values, second time: surfaceome filtered

jntvals1 =  read.table(file="clipboard", sep = "\t", header = T)

colnames(jntvals1)
head(jntvals1$JntComp1cp)
jntvalspos1 = jntvals1[jntvals1$JntComp1gene > 0 & jntvals1$JntComp1cp > 0,]
write.table(jntvalspos1,"jntcomp1_both_positive")
jntvalsneg1 = jntvals1[jntvals1$JntComp1gene < 0 & jntvals1$JntComp1cp < 0,]
write.table(jntvalsneg1,"jntcomp1_both_negative")

#jntcomp2

jntvals2 =  read.table(file="clipboard", sep = "\t", header = T)
jntvalspos2 = jntvals2[jntvals2$JntComp2gene > 0 & jntvals2$JntComp2cp > 0,]
write.table(jntvalspos2,"jntcomp2_both_positive")
jntvalsneg2 = jntvals2[jntvals2$JntComp2gene < 0 & jntvals2$JntComp2cp < 0,]
write.table(jntvalsneg2,"jntcomp2_both_negative")

#find the unique genes

jntvals1 =  read.table(file="clipboard", sep = "\t", header = T)
copyjntvals1 = jntvals1
rmd_jntvals1 = copyjntvals1[!duplicated(copyjntvals1$gene),]
write.table(rmd_jntvals1,"jnt1_transcriptional_targets_dupes_removed")

jntvals2 =  read.table(file="clipboard", sep = "\t", header = T)
copyjntvals2 = jntvals2
rmd_jntvals2 = copyjntvals2[!duplicated(copyjntvals2$gene),]
write.table(rmd_jntvals2,"jnt2_transcriptional_targets_dupes_removed")

#congyu's new data 
#129 of jnt2 and later 35 of jnt 1
neha = read.table(file="clipboard", sep = "\t", header = T) 
congyu = read.table(file="clipboard", sep = "\t", header = T)
intersection = inner_join(neha, congyu, by ="gene")
write.table(intersection, "substrates_only_FDR_from jnt1")

#check if the intersected data of jnt1 is included in the intersected data of jnt2
intersectedjnt1 = read.table(file ="clipboard", sep ="\t", header = T)
intersectedjnt2 = read.table(file ="clipboard", sep ="\t", header = T)
common_intersected_jnt12 = inner_join(intersectedjnt2, intersectedjnt1, by = "gene")

#Are the 35 also a subset of the 129?
jnt1_35 =  read.table(file="clipboard", sep = "\t", header = T) 
jnt2_129 =  read.table(file="clipboard", sep = "\t", header = T) 
common_intersected_jnt12 = inner_join(jnt2_129,jnt1_35 , by = "gene")
write.table(common_intersected_jnt12,"intersection_35jnt1_129jnt2")

#for 172 vs 172 of jnt comps 1 and 2
jnt1_172 =  read.table(file="clipboard", sep = "\t", header = T) 
jnt2_172 =  read.table(file="clipboard", sep = "\t", header = T) 
common_intersected_jnt12 = inner_join(jnt2_172,jnt1_172 , by = "gene")
write.table(common_intersected_jnt12,"intersection_172jnt1_172jnt2")

#for 272 vs congyu's substrates
all272 =  read.table(file="clipboard", sep = "\t", header = T)
congyu = read.table(file="clipboard", sep = "\t", header = T)
congyu272intersection = inner_join(all272, congyu, by = "gene")
write.table(congyu272intersection,"congyu_no_FD_intersect_all272")


# For the 172-129 data
leftover_jnt1 = read.table(file="clipboard", sep = "\t", header = T)
leftover_jnt2 = read.table(file="clipboard", sep = "\t", header = T)

colnames(leftover_jnt2)
#standardize the data

leftover_jnt2_std = leftover_jnt2 %>% mutate_at(c('JntComp2gene','JntComp2cp'), ~(scale(.) %>% as.vector))
leftover_jnt1_std = leftover_jnt1 %>% mutate_at(c('JntComp1gene','JntComp1cp'), ~(scale(.) %>% as.vector))

##for jnt2
#From the original transcriptome and proteome data, make a subset of the 43 genes and proteins
library(readxl)
mRNA_All_Results <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx", 
                               sheet = "Sheet4")

Protein_All_Results <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/Protein_All_Results.xlsx", 
                                  sheet = "Sheet5")

#find the subset of in leftover_jnt2
jnt2_leftovers = leftover_jnt2[,2]

rnaseq_leftover_jnt2 = mRNA_All_Results[mRNA_All_Results$gene %in% jnt2_leftovers,]
prot_leftover_jnt2 = Protein_All_Results[Protein_All_Results$gene %in% jnt2_leftovers,]
write.table(rnaseq_leftover_jnt2, "rnaseq_leftover_jnt2")
write.table(prot_leftover_jnt2,"prot_leftover_jnt2")
leftover_jnt2_original_data = inner_join(rnaseq_leftover_jnt2, prot_leftover_jnt2, by = "gene")
write.table(leftover_jnt2_original_data, "leftover_jnt2_original_data")

###############
#repeat for joint 1
jnt1_leftovers = leftover_jnt1[,2]

rnaseq_leftover_jnt1 = mRNA_All_Results[mRNA_All_Results$gene %in% jnt1_leftovers,]
prot_leftover_jnt1 = Protein_All_Results[Protein_All_Results$gene %in% jnt1_leftovers,]
write.table(rnaseq_leftover_jnt1, "rnaseq_leftover_jnt1")
write.table(prot_leftover_jnt1,"prot_leftover_jnt1")
leftover_jnt1_original_data = inner_join(rnaseq_leftover_jnt1, prot_leftover_jnt1, by = "gene")
write.table(leftover_jnt1_original_data, "leftover_jnt1_original_data")
################









#standardize these values

rnaseq_leftover_jnt2_std = rnaseq_leftover_jnt2[,-1] %>% mutate_all( ~(scale(.) %>% as.vector))
prot_leftover_jnt2_std = prot_leftover_jnt2[,-1:-2] %>% mutate_all( ~(scale(.) %>% as.vector))
#to innerjoin std values
rnalo2 = read.table(file="clipboard", sep = "\t", header = T)
protlo2 = read.table(file="clipboard", sep = "\t", header = T)
leftover_jnt2_original_data_std = inner_join(rnalo2, protlo2, by = "gene")
write.table(leftover_jnt2_original_data_std, "leftover_jnt2_original_data_std")























