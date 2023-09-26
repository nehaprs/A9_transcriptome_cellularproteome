#make subsets with only joint loading value 2 in all combinations
library(dplyr)
library(readr)


##secr_jnt_cp_jnt

#separate secretome values


#load the files with only the first joint loading values and all orthogonal loading values.

secr_jnt2_withCP <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/secr_cp/results/secr_jnt_cp_jnt/secr_joint_separated.csv", 
                      col_types = cols(V1 = col_skip()))

cp_joint2_withSECR <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/secr_cp/results/secr_jnt_cp_jnt/loading_values_yjoint.csv", 
                      col_types = cols(JntComp1_cp = col_skip()))


secr_orth <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/secr_cp/results/secr_jnt_cp_jnt/loading_values_xorth.csv")

cp_orth <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/secr_cp/results/secr_jnt_cp_jnt/loading_values_yorth.csv")
colnames(cp_orth)[1]<- "ID"

#inner_join

secr_jnt2_cp_jnt2 = inner_join(secr_jnt2_withCP,cp_joint2_withSECR,by="ID")
write.table(secr_jnt2_cp_jnt2,"secr_jnt2_cp_jnt2")

secr_jnt2_cp_orth = inner_join(secr_jnt2_withCP,cp_orth,by="ID")
write.table(secr_jnt2_cp_orth,"secr_jnt2_cp_orth")










##rnaseq_secretome
#load the files with only the first joint loading values and all orthogonal loading values.
rnaseq_jnt2_withSECR <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_sec/results/loading_values_xjoint.csv", 
                                 col_types = cols(V1 = col_skip()))
secr_jnt2_withRNA <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_sec/results/yjoint_separated.csv", 
                              col_types = cols(V1 = col_skip()))

rnaseq_orth_withSECR <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_sec/results/loading_values_xorth.csv")


secr_orth_withRNA <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_sec/results/loading_values_yorth.csv")

#inner_join
rnaseq_jnt2_secr_jnt2 = inner_join(rnaseq_jnt2_withSECR,secr_jnt2_withRNA,by='gene')
write.table(rnaseq_jnt2_secr_jnt2,"rnaseq_jnt2_secr_jnt2")
rnaseq_orth_secr_jnt2 = inner_join(rnaseq_orth_withSECR,secr_jnt2_withRNA,by='gene')
write.table(rnaseq_orth_secr_jnt2,"rnaseq_orth_secr_jnt2")





##rnaseq_cp
rnaseq_jnt2_withCP <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/loading_values_xjoint.csv", 
                               col_types = cols(V1 = col_skip()))
colnames(rnaseq_jnt2_withCP)[1] <-"gene"
rnaseq_orth_withCP <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/loading_values_xorth.csv")
cp_jnt2_withRNA <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/loading_values_yjoint.csv", 
                            col_types = cols(JntComp1_cp = col_skip()))
cp_orth_withRNA <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/loading_values_yorth.csv")
#inner_join
rnaseq_jnt2_cp_jnt2 = inner_join(rnaseq_jnt2_withCP,cp_jnt2_withRNA,by="gene")
write.table(rnaseq_jnt2_cp_jnt2,"rnaseq_jnt2_cp_jnt2")

rnaseq_orth_cp_jnt2 = inner_join(rnaseq_orth_withCP,cp_jnt2_withRNA,by="gene")
write.table(rnaseq_orth_cp_jnt2,"rnaseq_orth_cp_jnt2")
