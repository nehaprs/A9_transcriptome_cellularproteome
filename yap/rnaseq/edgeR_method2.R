#edgeR method2

library(edgeR)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

count_matrix <- as.matrix( read_excel("count_matrix.xlsx"))
counts = count_matrix[,-1]

#making sure counts has only numeric data
str(counts)
apply(counts, 2, function(x) all(is.numeric(x)))
counts <- apply(counts, 2, function(x) as.numeric(as.character(x)))
y	<- DGEList(counts=counts,	genes=count_matrix[,1])
head(y$genes)
rownames(y$counts)	<- rownames(y$genes)	<- y$genes$gene
y$genes$gene	<- NULL
y	<- calcNormFactors(y)
y$samples$group	=	factor(c("wt","wt","wt","yapko" ,"yapko" ,"yapko"))
plotMDS(y)
str(y$samples$group)
y	<- estimateDisp(y)
y$common.dispersion
plotBCV(y)

et	<- exactTest(y,	pair=c("wt", "yapko"))
summary(de	<- decideTestsDGE(et))
detags	<- rownames(y)[as.logical(de)]
plotSmear(et,	de.tags=detags) #this is the volcano plot showing DE genes
abline(h=c(-1,	1),	col="blue")

diffExpGenes	<- topTags(et,	n=10000,	p.value	=	0.05)
FCd = diffExpGenes$table[abs(diffExpGenes$table$logFC)>log2(1.5),]
write.table(FCd,	file="edgerDE_method2.txt",	sep	=	
              "\t",	row.names=TRUE,	col.names=NA)
yapde_fltrd <- read_excel("yapde_fltrd.xlsx")
comms = inner_join(yapde_fltrd, FCd, by = c("GeneName" = "genes"))
write_xlsx(comms,"common_deseq_edgeR2.xlsx")

edgerDE2 <- read_excel("BINF/multi-omics/yap/edgerDE2.xlsx")

library(readr)
a9_deg <- read_csv("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/DEG_filteredBy_FCnFDR.csv", 
                   col_types = cols(control1 = col_skip(), 
                                    control2 = col_skip(), control3 = col_skip(), 
                                    control4 = col_skip(), control5 = col_skip(), 
                                    control6 = col_skip(), ADAM9KD1 = col_skip(), 
                                    ADAM9KD2 = col_skip(), ADAM9KD3 = col_skip(), 
                                    ADAM9KD4 = col_skip(), ADAM9KD5 = col_skip(), 
                                    ADAM9KD6 = col_skip()))

names(edgerDE2)[names(edgerDE2) == 'genes'] <- 'GeneName'
a9_yap_common = inner_join(edgerDE2, a9_deg, by="GeneName")
library(writexl)
write_xlsx(a9_yap_common, "a9_yap_common_edgeR.xlsx")

a9_yapedg_drxn = data.frame()
a9_yap_common_edgeR <- read_excel("~/BINF/multi-omics/yap/a9_yap_common_edgeR.xlsx")

for(i in 1:nrow(a9_yap_common_edgeR)){
  row = a9_yap_common_edgeR[i,]
  if(sign(row$logFC.yap) == sign(row$logFC.a9)){
    a9_yapedg_drxn = rbind(a9_yapedg_drxn, row)
  }
}

write_xlsx(a9_yapedg_drxn,"a9_yap_edg_samedir.xlsx")
