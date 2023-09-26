#load log2FC of filtered data
mrna_log2= read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
prot_log2 = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)

library(OmicsPLS)
prot= impute_matrix(as.matrix(prot_log2))
mrna = impute_matrix(as.matrix(mrna_log2))

protT = t(prot) #this fn needs transpose of our original data frame
mrnaT = t(mrna)

#need to center data before cross validation, the tutorial asssumes all data centered around zero

mrnaTS = scale(mrnaT, center = TRUE, scale = TRUE)
mrnaTS = impute_matrix(mrnaTS)

#heatmap to make sure imputation of the data hasn't changed it much

#install.packages("gplots")
gplots::heatmap.2(cor(mrnaT,use = 'pair'), dendrogram='none', Rowv=F, Colv=F,trace='n',
                   col=gplots::bluered)
gplots::heatmap.2(cor(mrnaTS,use = 'pair'), dendrogram='none', Rowv=F, Colv=F,trace='n',
                  col=gplots::bluered)


protTS = scale(protT, center = TRUE, scale = TRUE)
protTS = impute_matrix(protTS)

crossval_o2m(mrnaTS,protTS,1:5,0:5,0:5,5)
#result: Minimal 5-CV error is at ax=1 ay=1 a=3
#result: Minimal 5-CV error is at ax=1 ay=0 a=3. Try this
#input_checker(mrnaTS,protTS) 
#results did not throw an error. So i/p ok
fit = o2m(mrnaTS,protTS,3,1,0)
summary(fit)
loadings = loadings(fit)

write.csv(loadings,"loading_values.csv")
scores = scores(fit)

write.csv(scores,"scores.csv")

#plot the loadings
dev.off()
plot(fit,"Xjoint",i=3,label = "colnames",size = 1.5)
plot(fit,"Yjoint", i=3,label = "colnames",size = 1.5)
plot(fit,"Xorth",label = "colnames",size = 1.5)
plot(fit,"Yorth",label = "colnames",size = 1.5)
plot(fit,"Xjoint", i = 1, j = 3, label = "colnames",size = 1.5)
plot(fit,"Xjoint", i = 1, j = 2, label = "colnames",size = 1.5)
plot(fit,"Xjoint", i = 2, j = 3, label = "colnames",size = 1.5)
plot(fit,"Yjoint", i = 1, j = 3, label = "colnames",size = 1)
plot(fit,"Yjoint", i = 1, j = 2, label = "colnames",size = 1.5)
plot(fit,"Yjoint", i = 2, j = 3, label = "colnames",size = 1.5)
#plot(Xjoint,Yjoint)
Xjoint
summary(fit)
plot(fit$Tt)
plot(fit$U)
yloadings = loadings(fit,"Yjoint")
write.csv(yloadings,"prot_joint_loading_values.csv")
scores(fit)

mrnaplot12 = plot(fit,"Xjoint", i = 1, j = 2, label = "colnames",size = 1.5)
protplot = plot(fit,"Yjoint", i = 1, j = 2, label = "colnames",size = 1.5)

# png("test.png")
# plot(fit,"Xjoint", i = 1, j = 2, label = "colnames",size = 1.5)
# par(new=TRUE)
# points("Yjoint", i = 1, j = 2)
# typeof(mrnaplot12)
# typeof(protplot)
# head(mrnaplot12)


#draw 3 dimensional plot of the joint loadings. Take csv of joint loadings, and 
#plot with x = v1, y = v2, z = v3
#did 3D plot in python. Much easier

#make heatmap of correlation ("correlogram")

install.packages("corrplot")
library(corrplot)

pearson_with_pvalue = read.table(file="clipboard", sep = "\t", header = T,row.names = 1)
