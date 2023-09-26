####lpa - edger
library(readxl)
edgerDE2 <- read_excel("BINF/multi-omics/yap/edgerDE2.xlsx")
protein_coding_genes <- read_excel("BINF/multi-omics/yap/protein_coding_genes.xlsx", 
                                     sheet = "Sheet1")

cds_in_edgerDE2 = edgerDE2[edgerDE2$genes %in% protein_coding_genes$gene,]
write.table(cds_in_edgerDE2, "cds_in_edgerDE2.xlsx")
library(dplyr)
library(readr)
a9_deg <- read_csv("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/DEG_filteredBy_FCnFDR.csv", 
                   col_types = cols(control1 = col_skip(), 
                                    control2 = col_skip(), control3 = col_skip(), 
                                    control4 = col_skip(), control5 = col_skip(), 
                                    control6 = col_skip(), ADAM9KD1 = col_skip(), 
                                    ADAM9KD2 = col_skip(), ADAM9KD3 = col_skip(), 
                                    ADAM9KD4 = col_skip(), ADAM9KD5 = col_skip(), 
                                    ADAM9KD6 = col_skip()))

common_cds_a9_edgeR_lpa = inner_join(a9_deg, cds_in_edgerDE2, by = c('GeneName' = 'genes') )

a9_yap_drxn = data.frame()
a9_yap_common = common_cds_a9_edgeR_lpa 

for(i in 1:nrow(a9_yap_common)){
  row = a9_yap_common[i,]
  if(sign(row$logFC.x) == sign(row$logFC.y)){
    a9_yap_drxn = rbind(a9_yap_drxn, row)
  }
}

#########
#starvation-edgeR

edger_starv_de <- read_excel("edger_starv_de.xlsx")

cds_in_edger_starv_de = edger_starv_de[edger_starv_de$genes %in% protein_coding_genes$gene,]

common_cds_a9_edgeR_starv = inner_join(a9_deg, cds_in_edger_starv_de, by = c('GeneName' = 'genes'))


a9_yap_drxn_st = data.frame()
a9_yap_common_st = common_cds_a9_edgeR_starv 

for(i in 1:nrow(a9_yap_common_st)){
  row = a9_yap_common_st[i,]
  if(sign(row$logFC.x) == sign(row$logFC.y)){
    a9_yap_drxn_st = rbind(a9_yap_drxn_st, row)
  }
}






write.table(a9_yap_drxn, 'a9_yap_drxn')
