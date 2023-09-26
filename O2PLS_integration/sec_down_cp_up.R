library(dplyr)

library(readr)
library(readxl)

#DEPs of secretomics filtered for FC and FDR
secretomics <- read_csv("sent/p_and_FC.csv")

#DEPs of cellular proteomics filtered for FC and FDR
cp <- read_csv("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/cell_proteomics_selected_for_FC_and_FDR.csv")

#finding common DEPs

sec_cp  = inner_join(secretomics,cp, by = "UniProt accession")
sec_cp <- sec_cp %>% rename_at("log2FC",~"log2FC.secretomics")
sec_cp <- sec_cp %>% rename_at("logFC",~"log2FC.Cell_Proteomics")
write.csv(sec_cp,"secretomics_cp_all_significant.csv")
#trying to check if the secretomics data has some issues

secr_pFC = secretomics
secr_dep = read.csv("clipboard")
ss = inner_join(secr_pFC, secr_dep, by = "UniProt accession")

secr_dep <- secr_dep %>% rename_at("UniProt.accession",~"UniProt accession")
#so 