setwd("~/BINF/multi-omics/yap/starvation")

#edgeR method2

library(edgeR)
library(dplyr)
library(AnnotationDbi)


count_matrix <- as.matrix(read_excel("count_matrix.xlsx"))
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

diffExpGenes	<- topTags(et,	n=20000,	p.value	=	0.05)
FCd = diffExpGenes$table[abs(diffExpGenes$table$logFC)>log2(1.5),]
write.table(FCd,	file="edgerDE_starv.txt",	sep	=	
              "\t",	row.names=TRUE,	col.names=NA)

library(readr)
a9_deg <- read_csv("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/DEG_filteredBy_FCnFDR.csv", 
                   col_types = cols(control1 = col_skip(), 
                                    control2 = col_skip(), control3 = col_skip(), 
                                    control4 = col_skip(), control5 = col_skip(), 
                                    control6 = col_skip(), ADAM9KD1 = col_skip(), 
                                    ADAM9KD2 = col_skip(), ADAM9KD3 = col_skip(), 
                                    ADAM9KD4 = col_skip(), ADAM9KD5 = col_skip(), 
                                    ADAM9KD6 = col_skip()))
edgerDE2 <- read_excel("edger_starv_de.xlsx")

names(edgerDE2)[names(edgerDE2) == 'genes'] <- 'GeneName'
a9_yap_common = inner_join(edgerDE2, a9_deg, by="GeneName")

library(writexl)
write_xlsx(a9_yap_common, "a9_yapstarv_common_edgeR.xlsx")

a9_yapedg_drxn = data.frame()
a9_yap_common_edgeR <- read_excel("a9_yap_starv_common.xlsx")

for(i in 1:nrow(a9_yap_common_edgeR)){
  row = a9_yap_common_edgeR[i,]
  if(sign(row$log2FC.yap) == sign(row$logFC.adam9)){
    a9_yapedg_drxn = rbind(a9_yapedg_drxn, row)
  }
}

nrow(a9_yap_common_edgeR)

write_xlsx(a9_yapedg_drxn,"a9_yap_satrv_edg_samedir.xlsx")
