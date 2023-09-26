library(pheatmap)
library(RColorBrewer)
#install.packages("DESeq2")
#install.packages("rlang")
library(rlang)
library(BiocManager)
library(DESeq2)
library(BiocManager)
#install.packages("tidyverse")
library(tidyverse)

library(readxl)

setwd("~/BINF/multi-omics/yap")

exp_counts <- read.table(file="clipboard", sep = "\t",
                         header = T, row.names = 1)
#matrix is count_matrix
exp_counts = exp_counts1
sample_data <- read_excel("metadata.xlsx")

data_deseq_group <- DESeqDataSetFromMatrix(countData =exp_counts ,
                                           colData = sample_data, design = ~ cell_type)

nrow(data_deseq_group)
data_deseq_group = data_deseq_group[rowSums(counts(data_deseq_group))>1]
#filtering the dataset for genes with 0 counts
nrow(data_deseq_group)

#from video
dds = DESeq(data_deseq_group)
res = results(dds)
BiocManager::install("apeglm")

library(apeglm)

reslfc = lfcShrink(dds, coef = 2)
(resOrdered <- res[order(res$padj),])
summary(res)

flt_vs_gc = as.data.frame(res$log2FoldChange)
head(flt_vs_gc)
write.csv(as.data.frame(resOrdered), file = "yapko_deseqDE.csv")

plotMA(reslfc, ylim = c(-1,1))

dds <- estimateSizeFactors(dds)
#se == std error
#se <- SummarizedExperiment(log2(counts(dds, normalize = TRUE)+ 1),colData = colData(dds))
se2 = rlog(dds, blind = FALSE)
se = se2
plotPCA(DESeqTransform(se),intgroup = "Sample")

#sample heatmap
sampleDists3 <- dist(t(assay(se)))
sampleDists3
sampleDistMatrix3 <- as.matrix(sampleDists3)
head(sampleDistMatrix3)
rownames(sampleDistMatrix3) <- paste( se$cell_type, se$Sample,
                                      sep="-" )

colnames(sampleDistMatrix3) <- paste( se$cell_type,se$Sample,
                                      sep="-" )
pheatmap(sampleDistMatrix3, clustering_distance_rows=sampleDists3,
         clustering_distance_cols=sampleDists3)

#filter 
yapko_deseqDE <- read_excel("yapko_deseqDE.xlsx")

yapde_sig = yapko_deseqDE[yapko_deseqDE$padj < 0.05,]

yapde_fltrd = yapde_sig[abs(yapde_sig$log2FoldChange) > log2(1.5),]

write_xlsx(yapde_fltrd, "yapde_fltrd.xlsx")

#compare it with our rnaseq deg

a9_deg <- read_csv("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/DEG_filteredBy_FCnFDR.csv", 
                                       col_types = cols(control1 = col_skip(), 
                                                                 control2 = col_skip(), control3 = col_skip(), 
                                                                 control4 = col_skip(), control5 = col_skip(), 
                                                                 control6 = col_skip(), ADAM9KD1 = col_skip(), 
                                                                  ADAM9KD2 = col_skip(), ADAM9KD3 = col_skip(), 
                                                                  ADAM9KD4 = col_skip(), ADAM9KD5 = col_skip(), 
                                                                  ADAM9KD6 = col_skip()))
names(yapde_fltrd)[names(yapde_fltrd) == '...1'] <- 'GeneName'

a9_yap = yapde_fltrd[yapde_fltrd$...1 %in% a9_deg$GeneName,]

a9_yap_common = inner_join(yapde_fltrd, a9_deg, by="GeneName")
library(writexl)
write_xlsx(a9_yap_common,"a9_yap_common.xlsx")

a9_yap_drxn = data.frame()


for(i in 1:nrow(a9_yap_common)){
  row = a9_yap_common[i,]
  if(sign(row$log2FoldChange) == sign(row$logFC)){
  a9_yap_drxn = rbind(a9_yap_drxn, row)
  }
}

write_xlsx(a9_yap_drxn,"a9_yap_samedir.xlsx")
