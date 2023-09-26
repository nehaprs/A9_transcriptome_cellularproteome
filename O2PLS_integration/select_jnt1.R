#make subsets with only joint loading value 1 in all combinations

library(dplyr)

##secr_jnt_cp_jnt
#load the files with only the first joint loading values and all orthogonal loading values.
secr_jnt1 <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/secr_cp/results/secr_jnt_cp_jnt/loading_values_xjoint.csv", 
                       col_types = cols(V2 = col_skip()))

cp_joint1 <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/secr_cp/results/secr_jnt_cp_jnt/loading_values_yjoint.csv", 
                    col_types = cols(JntComp2_cp = col_skip()))


secr_orth <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/secr_cp/results/secr_jnt_cp_jnt/loading_values_xorth.csv")

cp_orth <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/secr_cp/results/secr_jnt_cp_jnt/loading_values_yorth.csv")
colnames(cp_orth)[1]<- "ID"

#inner_join

secr_jnt1_cp_jnt1 = inner_join(secr_jnt1,cp_joint1,by="ID")
write.table(secr_jnt1_cp_jnt1,"secr_jnt1_cp_jnt1")

secr_jnt1_cp_orth = inner_join(secr_jnt1,cp_orth,by="ID")
write.table(secr_jnt1_cp_orth,"secr_jnt1_cp_orth")


##rnaseq_secretome
#load the files with only the first joint loading values and all orthogonal loading values.
rnaseq_jnt1_withSECR <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_sec/results/loading_values_xjoint.csv", 
                                 col_types = cols(V2 = col_skip()))
secr_jnt1_withRNA <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_sec/results/loading_values_yjoint.csv", 
                         col_types = cols(V2 = col_skip()))

rnaseq_orth_withSECR <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_sec/results/loading_values_xorth.csv")


secr_orth_withRNA <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_sec/results/loading_values_yorth.csv")

#inner_join
rnaseq_jnt1_secr_jnt1 = inner_join(rnaseq_jnt1_withSECR,secr_jnt1_withRNA,by='gene')
write.table(rnaseq_jnt1_secr_jnt1,"rnaseq_jnt1_secr_jnt1")
rnaseq_orth_secr_jnt1 = inner_join(rnaseq_orth_withSECR,secr_jnt1_withRNA,by='gene')
write.table(rnaseq_orth_secr_jnt1,"rnaseq_orth_secr_jnt1")

##rnaseq_cp
rnaseq_jnt1_withCP <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/loading_values_xjoint.csv", 
                                col_types = cols(V2 = col_skip()))
colnames(rnaseq_jnt1_withCP)[1] <-"gene"
rnaseq_orth_withCP <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/loading_values_xorth.csv")
cp_jnt1_withRNA <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/loading_values_yjoint.csv", 
                            col_types = cols(JntComp2_cp = col_skip()))
cp_orth_withRNA <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/loading_values_yorth.csv")
#inner_join
rnaseq_jnt1_cp_jnt1 = inner_join(rnaseq_jnt1_withCP,cp_jnt1_withRNA,by="gene")
write.table(rnaseq_jnt1_cp_jnt1,"rnaseq_jnt1_cp_jnt1")

rnaseq_orth_cp_jnt1 = inner_join(rnaseq_orth_withCP,cp_jnt1_withRNA,by="gene")
write.table(rnaseq_orth_cp_jnt1,"rnaseq_orth_cp_jnt1")
