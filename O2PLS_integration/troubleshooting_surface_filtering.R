library(dplyr)
#Find proteins in 91 but not in 86 

data86 =  read.table(file="clipboard", sep = "\t", header = T) 


in91_not_86 = anti_join(congyu91, data86)
in91_not_86 = in91_not_86[1:5, 1:5]


data_no_surf = read.table(file="clipboard", sep = "\t", header = T) 
surfaceome = read.table(file="clipboard", sep = "\t", header = T) 
surface_in_list = inner_join(data_no_surf, surfaceome)
surface_unique_in_list = surface_in_list %>% distinct()

#which are in 288 and not 272
in272 = read.table(file="clipboard", sep = "\t", header = T) 
in288_not_272 = anti_join(surface_unique_in_list, in272)
write.table(in288_not_272, "in288_not_272")
surf_not_fitered =  read.table(file="clipboard", sep = "\t", header = T) 
#find the 16 from this list
the15names = as.list(in288_not_272$gene)
components_of_the_15 = surf_not_fitered[surf_not_fitered$gene %in% the15names, ]
write.table(components_of_the_15, "components_of_the_15")

write.table(surface_in_list, "surface_in_list")
write.table(surface_unique_in_list, "surface_unique_in_list")

#of all 179 how many overlap with 91

all179 =  read.table(file="clipboard", sep = "\t", header = T)
all179in91 = inner_join(all179, congyu91)
write.table(all179in91, "all179in91")
