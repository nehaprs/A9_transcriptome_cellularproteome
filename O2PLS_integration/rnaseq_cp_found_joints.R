setwd("~/BINF/multi-omics")

#do o2pls without filtering for DEGs

library(readxl)
library(dplyr)
mRNA_All <- read_excel("ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx")

Protein_All <- read_excel("ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/Protein_All_Results.xlsx")


library(OmicsPLS)
#selecting for samples and genes alone

genes = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
#took only the first 5 samples of RNAseq as both genes and proteins need to have equal no. of samples
prots = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
#genes = mRNA_All[,c(1:13), sep = "\t", header = T,row.names = 1]. copy pasted since only that worked for scaling and centering for some reason.
#prots = Protein_All[,c(1,c(5:14))]
genes = impute_matrix(as.matrix(genes))

prots = impute_matrix(as.matrix(prots))
#samples as rows and genes/protein as columns
protT = t(prots)
geneT = t(as.matrix(genes))




#need to center data before cross validation, the tutorial asssumes all data centered around zero

geneTS = scale(geneT, center = TRUE, scale = TRUE)
geneTS = impute_matrix(geneTS)

protTS = scale(protT, center = TRUE, scale = TRUE)
protTS = impute_matrix(protTS)

                  
crossval_o2m_adjR2(geneTS,protTS,1:5,0:5,0:5,nr_folds = 5)
crossval_o2m(geneTS,protTS,1:3,0:2,0:2,5)
#Minimal 5-CV error is at ax=1 ay=1 a=2 
#Minimum MSE is 1.733012

fit = o2m(geneTS,protTS,2,1,1)
summary(fit)

loadings = loadings(fit)
loadings_xjoint = loadings(fit,"Xjoint")
write.csv(loadings_xjoint,"loading_values_xjoint.csv")
loadings_yjoint = loadings(fit,"Yjoint")
write.csv(loadings_yjoint,"loading_values_yjoint.csv")
loadings_xorth = loadings(fit,"Xorth")
write.csv(loadings_xorth,"loading_values_xorth.csv")
loadings_yorth = loadings(fit,"Yorth")
write.csv(loadings_yorth,"loading_values_yorth.csv")
loadings_gr_xjoint= loadings(fit,"gr_Xjoint")
write.csv(loadings_gr_xjoint,"loading_values_gr_xjoint.csv")

write.csv(loadings,"loading_values.csv")
scores = scores(fit)

write.csv(scores,"scores.csv")


#finding surface proteins 

#load proteins with joint components in rnaseq and cp

rnaseq_cp_all_joint <- read_csv("~/BINF/multi-omics/results for logFC/o2pls_all_data/rnaseq_cp/results/common_rnaseq_cp_all_joint.csv")
#secretome already loaded

rnaseq_cp_all_joint_surface = inner_join(rnaseq_cp_all_joint, surfaceome,by="ID")
write.csv(rnaseq_cp_all_joint_surface,"rnaseq_cp_all_joint_surface.csv")



#plotting joint components

plot(fit,"Xjoint",i=1,label = "colnames",size = 1.5)
#really big and messy plot of 1000s of values. Not worth plotting them all without filtering.