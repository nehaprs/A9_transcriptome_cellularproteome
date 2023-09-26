library(dplyr)
library(writexl)
congyutt = read.table(file="clipboard", sep = "\t", header = T) 
ttcomp1 =  read.table(file="clipboard", sep = "\t", header = T) 
ttcomp2 =  read.table(file="clipboard", sep = "\t", header = T)

comp1_in_congyutt = inner_join(congyutt, ttcomp1, by = "gene")
comp2_in_congyutt = inner_join(congyutt, ttcomp2, by = "gene")
write.table(comp1_in_congyutt, "comp1_in_congyu_tt")
write.table(comp2_in_congyutt, "comp2_in_congyu_tt")

overlap256_274 = inner_join(comp1_in_congyutt, comp2_in_congyutt, by = "gene")
write.table(overlap256_274, "overlap256_274")
not_in_overlap = anti_join(comp2_in_congyutt, overlap256_274, by = "gene")

overrrlap = comp1_in_congyutt[comp1_in_congyutt$gene %in% comp2_in_congyutt$gene &
                                comp1_in_congyutt$ID.x %in% comp2_in_congyutt$ID.x,]
write.table(overrrlap, "gene_and_ID_match256_vs_272")

overrrlap_gene = comp1_in_congyutt[comp1_in_congyutt$gene %in% comp2_in_congyutt$gene,]

overrrlap_geneRedo = inner_join(comp2_in_congyutt, comp1_in_congyutt, by = "gene")
write.table(overrrlap_geneRedo, "gene_only_match256v272")

#alternate way to find gene and ID match

gene_nID_256_272_redo = overrrlap_geneRedo[overrrlap_geneRedo$ID.x.x == overrrlap_geneRedo$ID.y.x &
                                             overrrlap_geneRedo$ID.y.x == overrrlap_geneRedo$ID.x.y &
                                             overrrlap_geneRedo$ID.x.y == overrrlap_geneRedo$ID.y.y,]

write.table(gene_nID_256_272_redo,"gene_nID_256_272_redo_unique_IDs")

###merge rows in ttcomp1 and ttcomp2 to find unique genes.
##writing function

merge_rows <- function(df){
  
  #create a dataframe to store the merged list
  merged_df <- data.frame(ID = character(), gene = character(), JntComp_gene = numeric(), JntComp_cp = numeric())
  
  #find unique values of genes 
  unique_values = unique(df$gene)
  
  #iterate over the unique values
  
  for(value in unique_values){
    rows <- df[df$gene == value,]
    
  #merge values in other columns, separated by comma
    merged_row <- data.frame(ID = paste(rows$ID, collapse = ", "),
                             gene = value,
                             JntComp_gene = paste(rows$JntComp2gene, collapse = ", "),
                             JntComp_cp = paste(rows$JntComp2cp, collapse = ", "))
    #append the merged_row to merged_df
    merged_df = rbind(merged_df, merged_row)
  }
  
  return(merged_df)
  
}

merged_ttcomp2 = merge_rows(ttcomp2)
write.table(merged_ttcomp2,"ttcomp2_repeatingGenesMerged")

common_mergedTT = inner_join(merged_ttcomp1, merged_ttcomp2, by = "gene")
write.table(common_mergedTT,"mergedTT_commonJn1jn2")

##addition after completing PTTs
#finding union of jn1 and jn2 candidate tts

tt_jnt1 =  read.table(file="clipboard", sep = "\t", header = T)
tt_jnt2 =  read.table(file="clipboard", sep = "\t", header = T)
tt_jnt1n = tt_jnt1[,c(2,3)]
tt_jnt2n = tt_jnt2[,c(2,3)]

tt_union = full_join(tt_jnt1n ,tt_jnt2n )
write_xlsx(tt_union,"tt_jnts_union.xlsx")
#remove repeated peptides from tt_union
tt_un_repeats = subset(tt_union, grepl(",", ID))
tt_un_repeats_ID = tt_un_repeats[,1]
write.table(tt_un_repeats_ID, "clipboard", sep="\t", row.names=FALSE)
###from these, find the unreviewed entries manually 
#upload the unreviewed list
unreviewed_ids <- read_excel("idmapping_reviewed_false_2023_07_21 (1).xlsx")
#remove them from tt_union
tt_union2 = data.frame(tt_union, stringsAsFactors = FALSE)
unrvd2 = data.frame(unreviewed_ids, stringsAsFactors = FALSE)

tt_jnts_union <- read_excel("tt_jnts_union.xlsx", 
                               sheet = "Sheet4")
#from tt_union, how many proteins are significantly DE'd on A9KD
#load FC data
Protein_All_Results <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/Protein_All_Results.xlsx")
tt_union_proteome = Protein_All_Results[Protein_All_Results$ID %in% tt_jnts_union$ID, c(1,3,15,21)]
nn = anti_join( tt_union_proteome,tt_jnts_union)
tt_union_proteomeSig = tt_union_proteome[tt_union_proteome$adj.P.Val < 0.05,]
tt_neg = tt_union_proteomeSig[tt_union_proteomeSig$logFC < 0,]
tt_pos = tt_union_proteomeSig[tt_union_proteomeSig$logFC > 0,]

write_xlsx(tt_neg,"tt_neg.xlsx")
write_xlsx(tt_pos,"tt_pos.xlsx")


#different fold change cut-offs for TTs

#neg_TTs
tt_neg1.2 = tt_neg[tt_neg$logFC < -1.2,]
write_xlsx(tt_neg1.2,"tt_neg1.2.xlsx")

tt_neg1 = tt_neg[tt_neg$logFC < -1,]
write_xlsx(tt_neg1,"tt_neg1.xlsx")

tt_neg1.5 = tt_neg[tt_neg$logFC < -1.5,]
write_xlsx(tt_neg1.5,"tt_neg1.5.xlsx")

tt_neg0.5 = tt_neg[tt_neg$logFC < -0.5,]
write_xlsx(tt_neg0.5,"tt_neg0.5.xlsx")

#pos_TTs
tt_pos1.2 = tt_pos[tt_pos$logFC > 1.2,]
write_xlsx(tt_pos1.2,"tt_pos1.2.xlsx")

tt_pos1 = tt_pos[tt_pos$logFC > 1,]
write_xlsx(tt_pos1,"tt_pos1.xlsx")

tt_pos1.5 = tt_pos[tt_pos$logFC > 1.5,]
write_xlsx(tt_pos1.5,"tt_pos1.5.xlsx")

tt_pos0.5 = tt_pos[tt_pos$logFC > 0.5,]
write_xlsx(tt_pos0.5,"tt_pos0.5.xlsx")
