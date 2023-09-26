##finding surface proteins
library(readxl)
library(dplyr)
surfaceome <- read_excel("~/BINF/secretomics/references/table_S3_surfaceome.xlsx", 
                                     sheet = "Sheet1")

dataaa = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)


comb = inner_join(dataaa, surfaceome,by="ID")
write.table(comb,"rnaseq_orth_withCP_surface")

##rnaseq orth-cp up

rnaseq = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)

cp =  read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
joint_rnaseq_cp = inner_join(rnaseq,cp, by = "gene")
write.table(joint_rnaseq_cp,"rnaseq_orth_cp_jnt")

##select large joint values

dataaa = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)

cp_positives = subset(dataaa, dataaa$JntComp2_cp > 0)
install.packages("formattable")
library(formattable)
#if jntComp_cp > 0.01, color blue, else grey
cp_positives_large = formattable(cp_positives,list(JntComp2_cp = formatter("span",style = ~style(display = "block",
                                                                       font.weight = "bold", 
                                                                       color = "white",
                                                                       "border-radius" = "4px",
                                                                       "padding-right" = "4px",
                                                                       "background-color" =  
                                                                         ifelse(JntComp2_cp > 0.01,"darkblue",
                                                                          ifelse(JntComp2_cp < 0.01,"grey",
                                                                                 ifelse(JntComp2_cp == 0.01,"grey",NA)))))))

both_large = formattable(cp_positives,list(rnaseq_orth = formatter("span",style = ~style(display = "block",
                                                                            font.weight = "bold", 
                                                                            color = "white",
                                                                            "border-radius" = "4px",
                                                                            "padding-right" = "4px",
                                                                            "background-color" =  
                                                                              ifelse(JntComp2_cp > 0.01 & abs(rnaseq_orth)>0.01,"darkblue",
                                                                                ifelse("grey",NA))))))
write.table(both_large,"both_large")

cp_positives = cp_positives %>% arrange(rnaseq_orth,JntComp2_cp)
write.table(cp_positives,"rnaseq_orth_cp_posJnt_ordered")
