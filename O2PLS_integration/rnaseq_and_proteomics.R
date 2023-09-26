
# select FC > 1.5, p-value <0.05

mrna_all = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
prot_all = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
mrna_FC = mrna_all[abs(mrna_all$logFC) > abs(log2(1.5)),]
mrna_FC_FDR = mrna_FC[abs(mrna_FC$FDR)< 0.05,]
write.csv(mrna_FC_p,'rnaseq_selected_for_FC_and_p.csv')

prot_all = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
prot_FC = prot_all[abs(prot_all$logFC) > abs(log2(1.5)),]
prot_FC_FDR = prot_FC[abs(prot_FC$adj.P.Val)< 0.05,]
write.csv(prot_FC_FDR,'cell_proteomics_selected_for_FC_and_FDR.csv')


#from there written csvs, columns containing names and replicates si1-si5 and kd1- kd5
#are used for o2pls analysis 

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_201')
library("rJava")
library("Deducer")


library(OmicsPLS)
install.packages('xlsx')     
library(xlsx)   
install.packages('readxl')     
library(readxl)
install.packages('writexl')
library(writexl)


#install.packages("OmicsPLS")
#install.packages("parallel")
#install.packages("ggplot2")

mrna = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
prot = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
secr = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)


prot= impute_matrix(as.matrix(prot))
protT = t(prot) #this fn needs transpose of our original data frame

mrna  = impute_matrix(as.matrix(mrna))
mrnaT = t(mrna)

secr  = impute_matrix(as.matrix(secr))
secrT = t(secr)

#need to center data before cross validation, the tutorial asssumes all data centered around zero

mrnaTS = scale(mrnaT, center = TRUE, scale = TRUE)
mrnaTS = impute_matrix(mrnaTS)

protTS = scale(protT, center = TRUE, scale = TRUE)
protTS = impute_matrix(protTS)

secrTS = scale(secrT, center = TRUE, scale = TRUE)
secrTS = impute_matrix(secrTS)

mrnaTS_df = as.data.frame(mrnaTS)
protTS_df = as.data.frame(protTS)
secrTS_df = as.data.frame(secrTS)


write_xlsx(mrnaTS_df,"mrnaTS.xlsx")
write_xlsx(protTS_df,"protTS.xlsx")
write_xlsx(secrTS_df,"secrTS.xlsx")


write.table(mrnaTS,"mrnaTS")
write.table(protTS,"protTS")
write.table(secrTS,"secrTS")
#heatmap to make sure imputation of the data hasn't changed it much

install.packages("gplots")
gplots::heatmap.2(cor(protT,use = 'pair'), dendrogram='none', Rowv=F, Colv=F,trace='n',
                  breaks=seq(-1,1,length.out = 25), col=gplots::bluered)
gplots::heatmap.2(cor(protTS,use = 'pair'), dendrogram='none', Rowv=F, Colv=F,trace='n',
                  breaks=seq(-1,1,length.out = 25), col=gplots::bluered)
#both heatmaps look very similar, so imputation hasn't made a lot of difference
#only prot data imputed. there was no missing data in rnaseq


crossval_o2m_adjR2(mrnaTS,protTS,1:6,1:6,1:6,2)
# result: min at n =3, nx = ny =1
# debugging: find if you need to increase the upper limit of the 3 vectors used
fit = o2m(mrnaTS,protTS,3,1,1)
summary(fit)
plot(fit,"Xjoint",label = "colnames",size = 1.5)
plot(fit,"Yjoint",label = "colnames",size = 1.5)
plot(fit,"Xorth",label = "colnames",size = 1.5)
plot(fit,"Xorth",label = "colnames",size = 1.5)
plot(fit,"Yjoint", i = 1, j = 2, label = "colnames",size = 1.5)
#plotting of correlation

#install.packages("magrittr")
#install.packages("ggplot2")
#install.packages("gridExtra")
#nstall.packages("stringr")
#install.packages("gplots")
#install.packages("reshape2")

library(magrittr)
library(ggplot2)
library(gridExtra)
BiocManager::install("illuminaHumanv3.db")
library(illuminaHumanv3.db)



# Color names
LLmodule <- c("ILMN_1690209",'ILMN_1766551', 'ILMN_1749131', 'ILMN_1688423',
              'ILMN_2102670', 'ILMN_1792323', 'ILMN_1899034', 'ILMN_1806721',
              'ILMN_1695530', 'ILMN_1726114', 'ILMN_1751625', 'ILMN_1726114',
              'ILMN_1753648', 'ILMN_1779043')
LLnr <- which(colnames(mrnaTS) %in% LLmodule)
#rna_genenames <- select(illuminaHumanv3.db,
                        #keys = colnames(rna)[LLnr],
                        #keytype = "PROBEID", columns = "SYMBOL")[,2]

#Plot loadings with OmicsPLS plot method ###
p_prot <- plot(fit, loading_name="Yj", i=1, j=2, label="c", alpha=1)
p_prot


head(prot[,(5+5)])

##repeating the above on log2FC



# r file rnaseq_prot_logFC
