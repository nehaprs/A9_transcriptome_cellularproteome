library(dplyr)

#load FC data of the leftover 43 genes


fcdata = read.table(file="clipboard", sep = "\t", header = T) 

#make df of rna and prot data

fcrna = fcdata[,c(1:9)]
fcprot = fcdata[,c(1:3, 10:14)]



fcrna_std = fcrna[,-c(1:3)] %>% mutate_all(~(scale(.) %>% as.vector))
fcprot_std = fcprot[,-c(1:3)] %>% mutate_all(~(scale(.) %>% as.vector))
#extract names
fcrna_names = fcrna[,c(1,2,3)]
fcprot_names = fcprot[,c(1,2,3)]
#combine names and standardized values
fcrna_std_named = cbind(fcrna_names, fcrna_std)
fcprot_std_named = cbind(fcprot_names, fcprot_std)
#combine the standardized values to form 1 dataframe with standardized values
library(dplyr)
fcdata_std = inner_join(fcrna_std_named, fcprot_std_named, by ="gene")

# declaring an empty data frame to populate gene name, pearson, p-value
leftover43_jnt2_pearson = data.frame(gene = character(), PearsonCoeff = numeric(), p_value = numeric())

#find pearson's correlation gene by gene
head(fcdata_std[,3])
for(i in 1: nrow(fcdata_std))
{rowi = fcdata_std[i,]
leftover43_jnt2_pearson[i,1] = rowi[3]
fcrnai = unlist(as.vector(rowi[,4:8]))
fcproti = unlist( as.vector(rowi[,12:16]))
pearsi = cor.test(fcrnai, fcproti, method = "pearson")
leftover43_jnt2_pearson[i,2] = pearsi$estimate
leftover43_jnt2_pearson[i,3] = pearsi$p.value

} 


write.table(leftover43_jnt2_pearson, "leftover43_jnt2_pearson")

###################################

##same things for joint component 1
#load FC data of the leftover 43 genes

fcdata = read.table(file="clipboard", sep = "\t", header = T) 

#make df of rna and prot data

fcrna = fcdata[,c(1:9)]
fcprot = fcdata[,c(1:3, 10:14)]



fcrna_std = fcrna[,-c(1:3)] %>% mutate_all(~(scale(.) %>% as.vector))
fcprot_std = fcprot[,-c(1:3)] %>% mutate_all(~(scale(.) %>% as.vector))
#extract names
fcrna_names = fcrna[,c(1,2,3)]
fcprot_names = fcprot[,c(1,2,3)]
#combine names and standardized values
fcrna_std_named = cbind(fcrna_names, fcrna_std)
fcprot_std_named = cbind(fcprot_names, fcprot_std)
#combine the standardized values to form 1 dataframe with standardized values
library(dplyr)
fcdata_std = inner_join(fcrna_std_named, fcprot_std_named, by ="gene")

# declaring an empty data frame to populate gene name, pearson, p-value
leftover43_jnt1_pearson = data.frame(gene = character(), PearsonCoeff = numeric(), p_value = numeric())

#find pearson's correlation gene by gene

for(i in 1: nrow(fcdata_std))
{rowi = fcdata_std[i,]
leftover43_jnt1_pearson[i,1] = rowi[3]
fcrnai = unlist(as.vector(rowi[,4:8]))
fcproti = unlist( as.vector(rowi[,12:16]))
pearsi = cor.test(fcrnai, fcproti, method = "pearson")
leftover43_jnt1_pearson[i,2] = pearsi$estimate
leftover43_jnt1_pearson[i,3] = pearsi$p.value

} 


write.table(leftover43_jnt1_pearson, "leftover43_jnt1_pearson")


##########################
#Pearson against adam9 transcriptomics

#get tanscriptome fc of the 43

fc_43_trans = read.table(file="clipboard", sep = "\t", header = T)
fc_43_trans = fc_43_trans[,1:8]
#scale
fc_43_trans_std = fc_43_trans[,-c(1:2)] %>% mutate_all(~(scale(.) %>% as.vector))
fc_43_trans_stdNamed[1,]
fc43names = fc_43_trans[,c(1,2)]
fc_43_trans_stdNamed = cbind(fc43names,fc_43_trans_std)
#declare a df to write pearsons
# declaring an empty data frame to populate gene name, pearson, p-value
leftover43_jnt1A9_pearson = data.frame(gene = character(), PearsonCoeff = numeric(), p_value = numeric())
#pearson
for(i in 2: nrow(fc_43_trans_stdNamed))
{ a9 = fc_43_trans_stdNamed[1,]
  rowi = fc_43_trans_stdNamed[i,]
leftover43_jnt1A9_pearson[i,1] = rowi[2]
fcrnai = unlist(as.vector(rowi[,-c(1,2)]))
a9fc = unlist( as.vector(a9[,-c(1,2)]))
pearsi = cor.test(fcrnai, a9fc, method = "pearson")
leftover43_jnt1A9_pearson[i,2] = pearsi$estimate
leftover43_jnt1A9_pearson[i,3] = pearsi$p.value

} 

write.table(leftover43_jnt1A9_pearson,"pearson_jn1_A9_transcript")


#combined data of the earlier 129 and 37 rescued by pearson


data167 = read.table(file="clipboard", sep = "\t", header = T) 
congyu91 = read.table(file="clipboard", sep = "\t", header = T)
intersect_167_subs_91_subs = inner_join(data167, congyu91)

write.table(intersect_167_subs_91_subs, "intersect_167_subs_91_subs")


