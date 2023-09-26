setwd("~/BINF/multi-omics")
#First part: finding common genes/proteins without O2PLS. Second part uses O2PLS
library(readxl)
all_genes <- read_excel("ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx")
DEG = all_genes[abs(all_genes$logFC) > log2(1.5),]
DEG = DEG[DEG$FDR < 0.05,]
#write.csv(DEG,"DEG_filteredBy_FCnFDR.csv")

#finding all DEG that doesn't have any significant change

no_change_genes = subset(all_genes, !(all_genes$GeneName %in% DEG$GeneName))
write.csv(no_change_genes,"notDEGs.csv")

#upload secretome data

library(readr)
secretome <- read_csv("~/BINF/secretomics/sent/p_and_FC.csv")
secretome_down = secretome[secretome$log2FC < 0,]
#write.csv(secretome_down,"downregulated_secretome.csv")

#find commons between not DE genes and downregulated secretome

sec_down_transc_noDEG  = subset(no_change_genes, no_change_genes$GeneName %in% secretome_down$gene_name)


#include gene names to secretome dataset
#done in python, colab file multi-omics under the text documentation saying this 

#write.csv(sec_down_transc_noDEG,"sec_down_rnaseq_noDE.csv")

#finding surface proteins 
library(readxl)
surfy <- read_excel("~/BINF/secretomics/references/table_S3_surfaceome.xlsx", 
                     +     sheet = "in silico surfaceome only")
library(dplyr)
surfy <- surfy %>% rename_at("UniProt gene",~"UniProt_gene")

in_surfy = subset(sec_down_transc_noDEG, sec_down_transc_noDEG$GeneName %in% surfy$UniProt_gene)

 #write.csv(in_surfy,"sec_down_rnaseq_noDE_surfy.csv")
 
#x = genes
#y = secretome

library(OmicsPLS)
library(readxl)
secretome_all <- read_excel("~/BINF/secretomics/intermediate_outputs/secretome_all.xlsx")
secretome_all= impute_matrix(as.matrix(secretome_all))
genes_all = impute_matrix(as.matrix(all_genes))

secT = t(secretome_all)
geneT = t(genes_all)

#need to center data before cross validation, the tutorial asssumes all data centered around zero

secTS = scale(secT, center = TRUE, scale = TRUE)
mrnaTS = impute_matrix(mrnaTS)
