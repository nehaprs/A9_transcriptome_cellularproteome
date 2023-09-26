library(readr)
all_corrs <- read_csv("BINF/multi-omics/results for logFC/excels/pearson_filtered_significant_both_directions.csv")

prot_names2 = sub( '_HUMAN'," ", all_corrs$Protein)
head(prot_names2)
all_corrs$prot_names = prot_names2

chh = c(1,2,3,4)
chh[1]
for (i in (1: nrow(all_corrs))){
 protein_ith = all_corrs[i,6]
 rna_ith = all_corrs[i,2]
 if(protein_ith == rna_ith){
   print(all_corrs[i,])
 }
 
}

# result: not a single rna-protein pair has a filtered significant interaction.
#So we need to try with the full unfiltered list
corrs_no_filter <- read_csv("BINF/multi-omics/results for logFC/excels/pearsons_correlation_rna_prot.csv")

prot_names = sub( '_HUMAN',"", corrs_no_filter$Protein)
corrs_no_filter$prot_names = prot_names

comms = list()
for (i in (1: nrow(corrs_no_filter))){
  protein_ith = corrs_no_filter$prot_names[i]
  #print(protein_ith)
  rna_ith = corrs_no_filter$RNA[i]
  #print(rna_ith)
  if (is.na(rna_ith)){
    print('Missing')
  }
  else if(protein_ith == rna_ith){
    comms = append(comms,corrs_no_filter[i,])
  
  }
  
}

write.csv(comms,"common_prot_rna.csv")

library(readxl)

gene_names_of_mapped_prots <- read_excel("BINF/multi-omics/results for logFC/excels/gene_names_of_mapped_prots.xlsx")

gene_name = gene_names_of_mapped_prots$`Gene Names`

if ('ADAM9'  %in% corrs_no_filter$RNA){print("y")}

num = 0
genes_same_as_encoders = list()
#gene_names_list <- read_csv("BINF/multi-omics/results for logFC/excels/list_gene_names.csv")
for (j in (1:nrow(gene_names_of_mapped_prots))){
for (i in (3:ncol(gene_names_of_mapped_prots))){
 
  if (gene_names_of_mapped_prots[j,i] %in% corrs_no_filter$RNA ){
    print(gene_names_of_mapped_prots[j,i])
    genes_same_as_encoders = append(genes_same_as_encoders,gene_names_of_mapped_prots[j,i])
  }
  
}
}
num

length(corrs_no_filter$RNA)
