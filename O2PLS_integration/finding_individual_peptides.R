library(dplyr)

#load comma separated files for tt1 and tt2 w
#we have them in the workspace as merged_ttcomp1 and mer

#repeated genes

id_repeated = c('P17301', 'E7ESP4', 'Q96EB3', 'P68104', 'Q9UJU6', 'B4DDD6', 'Q6IPJ9', 'O00515')

peptides_cellular_proteomics <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/peptides_cellular_proteomics.xlsx")

#create a new df for individual peptides

individual_peptides = data.frame(ID = character(), peptide = character())

#get the peptides from the datasets

for(i in 1:8){
  id_i = id_repeated[i]
  rows = peptides_cellular_proteomics[peptides_cellular_proteomics$ID == id_i,]
  individual_peptides = rbind(individual_peptides, rows)
}



write.table(individual_peptides,"individual_peptides_256n272")

#for all tt of jntComp1

ids_jnt1TT <- data.frame(Column1 = unlist(strsplit(as.character(merged_ttcomp1$ID), ", ")))

ids_jnt1TT = merged_ttcomp1[grepl(",", merged_ttcomp1$ID),]

##unstring ids in ids_jnt1TT

#create an empty df

unstr_1TT = data.frame(ID = character(), gene = character(), stringsAsFactors = FALSE)

#iterate

for(i in 1:nrow(ids_jnt1TT)){
  id_list = ids_jnt1TT$ID
  xx = strsplit(ids_jnt1TT$ID[i], ", ")[[1]]
  for(j in 1:length(xx)){ 
    # Create new rows in df2
    xxj = xx[j]
    new_rows <- data.frame(ID = xxj,
                           gene = rep(ids_jnt1TT$gene[i], length(xxj)),
                           stringsAsFactors = FALSE)
    
    unstr_1TT <- rbind(unstr_1TT, new_rows)}
 
}

write.table(unstr_1TT,"list_of_repeating_IDs_allTTJn1")

##################

#find individual peptides for all jnt1 IDs



id_repeated = unstr_1TT$ID

peptides_cellular_proteomics <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/peptides_cellular_proteomics.xlsx")

#create a new df for individual peptides

individual_peptides = data.frame(ID = character(), peptide = character())

#get the peptides from the datasets

for(i in 1:length(id_repeated)){
  id_i = id_repeated[i]
  rows = peptides_cellular_proteomics[peptides_cellular_proteomics$ID == id_i,]
  
  individual_peptides = rbind(individual_peptides, rows)
}

which(id_repeated == 'P08133', arr.ind = T)



##add the gene column to the df individual peptides

gene_col = c()

#populate gene by gene names in unstr_id

for(i in 1: length(individual_peptides$ID)){
  ID_name = individual_peptides$ID[i]
  for(j in 1: nrow(unstr_1TT)){
    if(unstr_1TT$ID[j]== ID_name){
      gene_name = unstr_1TT$gene[j]
      gene_col = append(gene_col, gene_name)
    }
  }
  
}

individual_peptides = cbind(individual_peptides, gene_col)

write.table(individual_peptides,"individual_peptides_allJnt1")


#######
##individual peptides for joint component 2
ttcomp2 =  read.table(file="clipboard", sep = "\t", header = T)
merged_ttcomp2 <- read_excel("ttcomp2_repeatingGenesMerged.xlsx")
ids_jnt2TT <- data.frame(Column1 = unlist(strsplit(as.character(merged_ttcomp2$ID), ", ")))

ids_jnt2TT = merged_ttcomp2[grepl(",", merged_ttcomp2$ID),]

##unstring ids in ids_jnt1TT

#create an empty df

unstr_2TT = data.frame(ID = character(), gene = character(), stringsAsFactors = FALSE)

#iterate

for(i in 1:nrow(ids_jnt2TT)){
  id_list = ids_jnt2TT$ID
  xx = strsplit(ids_jnt2TT$ID[i], ", ")[[1]]
  for(j in 1:length(xx)){ 
    # Create new rows in df2
    xxj = xx[j]
    new_rows <- data.frame(ID = xxj,
                           gene = rep(ids_jnt2TT$gene[i], length(xxj)),
                           stringsAsFactors = FALSE)
    
    unstr_2TT <- rbind(unstr_2TT, new_rows)}
  
}

write.table(unstr_2TT,"list_of_repeating_IDs_allTTJn2")

##################

#find individual peptides for all jnt2 IDs



id_repeated = unstr_2TT$ID

peptides_cellular_proteomics <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEPs/peptides_cellular_proteomics.xlsx")

#create a new df for individual peptides

individual_peptides = data.frame(ID = character(), peptide = character())

#get the peptides from the datasets

for(i in 1:length(id_repeated)){
  id_i = id_repeated[i]
  rows = peptides_cellular_proteomics[peptides_cellular_proteomics$ID == id_i,]
  
  individual_peptides = rbind(individual_peptides, rows)
}




##add the gene column to the df individual peptides

gene_col = c()

#populate gene by gene names in unstr_id

for(i in 1: length(individual_peptides$ID)){
  ID_name = individual_peptides$ID[i]
  for(j in 1: nrow(unstr_2TT)){
    if(unstr_2TT$ID[j]== ID_name){
      gene_name = unstr_2TT$gene[j]
      gene_col = append(gene_col, gene_name)
    }
  }
  
}

individual_peptides = cbind(individual_peptides, gene_col)

write.table(individual_peptides,"individual_peptides_allJnt2")
