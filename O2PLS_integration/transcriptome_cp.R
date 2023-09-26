#of the transcriptome-cellular proteome joint components, find the prots with the highest value for each joint component

sample_data = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
#sample_data + secr_cp_common_joint_components
library(dplyr)
library(tidyr)

#rnaseq joint component 1
scr1 = sample_data[,c(1,2)]
scr1_srt = scr1[sort(abs(scr1$"V1gene"),decreasing=T,index.return=T)[[2]],]
scr1_srt_top100 = head(scr1_srt,100)
scr1_srt_top100_ordered = scr1_srt_top100[sort((scr1_srt_top100$"V1gene"),decreasing=F,index.return=T)[[2]],]
write.csv(scr1_srt_top100_ordered,"rnaseq1_srt_top100_ordered.csv")

#rnaseq joint component 2
scr2 = sample_data[,c(1,3)]
scr2_srt = scr2[sort(abs(scr2$"V2gene"),decreasing=T,index.return=T)[[2]],]
scr2_srt_top100 = head(scr2_srt,100)
scr2_srt_top100_ordered = scr2_srt_top100[sort((scr2_srt_top100$"V2gene"),decreasing=F,index.return=T)[[2]],]
write.csv(scr2_srt_top100_ordered,"scr2_srt_top100_ordered.csv")

#cp joint component 1
cp1 = sample_data[,c(1,4,5)]
cp1_srt = cp1[sort(abs(cp1$"V1cp"),decreasing=T,index.return=T)[[2]],]
cp1_srt_top100 = head(cp1_srt,100)
cp1_srt_top100_ordered = cp1_srt_top100[sort((cp1_srt_top100$"V1cp"),decreasing=F,index.return=T)[[2]],]
write.csv(cp1_srt_top100_ordered,"cp1_srt_top100_ordered.csv")

#cp joint component 2
cp2 = sample_data[,c(1,4,6)]
cp2_srt = cp2[sort(abs(cp2$"V2cp"),decreasing=T,index.return=T)[[2]],]
cp2_srt_top100 = head(cp2_srt,100)
cp2_srt_top100_ordered = cp2_srt_top100[sort((cp2_srt_top100$"V2cp"),decreasing=F,index.return=T)[[2]],]
write.csv(cp2_srt_top100_ordered,"cp2_srt_top100_ordered.csv")