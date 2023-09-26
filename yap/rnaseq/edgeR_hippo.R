library(BiocManager)

#BiocManager::install("edgeR")
library(edgeR)
library(readxl)
count_matrix <- as.matrix( read_excel("count_matrix.xlsx"))
DataGroups <- c("wt","wt","wt","yapko" ,"yapko" ,"yapko")
d <- DGEList(counts=count_matrix,group=factor(DataGroups))
#remove genenames from count data
counts = count_matrix[,-1]
head(counts)
#making sure counts has only numeric data
str(counts)
apply(counts, 2, function(x) all(is.numeric(x)))
counts <- apply(counts, 2, function(x) as.numeric(as.character(x)))

d <- DGEList(counts=counts,group=factor(DataGroups))
d
dim(d)
d.full <- d # keep the old one in case we mess up
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum) # total gene counts per sample
#kept only those with cpm abouve 0.01
keep <- rowSums(cpm(d)>0.01) >= 2
d <- d[keep,]
dim(d)

d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d)
d

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

d1 <- estimateCommonDisp(d, verbose=T)

names(d1)

d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

et12 <- exactTest(d1, pair=c(1,2))
topTags(et12, n=10)
diffExpGenes	<- topTags(et12,	n=10000,	p.value	=	0.05)
de1 <- decideTestsDGE(et12, adjust.method="BH")
write.table(diffExpGenes$table,	file="yap_edger.txt",	sep	=	
              "\t",	row.names=TRUE,	col.names=NA)
# differentially expressed tags from the naive method in d1
de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")

counts_numbered <- read_excel("counts_numbered.xlsx")
yap_edger <- read_excel("yap_edger.xlsx")
yap_edger_final = inner_join(yap_edger, counts_numbered, by ="...1")
write_xlsx(yap_edger_final,"yap_edger_named.xlsx")

#compare with deseq
yapde_fltrd <- read_excel("yapde_fltrd.xlsx")

commons = yapde_fltrd[yapde_fltrd$GeneName %in% yap_edger_final$Column1,]
write_xlsx(commons,"common_edger_deseq.xlsx")
