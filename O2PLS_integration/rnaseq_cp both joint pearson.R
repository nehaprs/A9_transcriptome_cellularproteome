############
#midnight: rnaseq_cp both joint pearson
############

library(readr)
gene <- read_csv("~/BINF/multi-omics/o2pls_all_data/rnaseq_cp/results/loading_values_xjoint.csv")




cprot <- read_csv("~/BINF/multi-omics/o2pls_all_data/rnaseq_cp/results/y_joint_loading_values_for_python.csv")

gene2 <- data.frame(name = gene$...1,
                    gene_loading = gene$V1)

###create 2d vectors from loading values

#gene2$gene_loading <- as.vector(gene2$gene_loading)



cprot2 <- data.frame(name = cprot$...1, cprot_loading = cprot$JntComp1_cp)
#cprot2$cprot_loading <- as.vector(cprot2$cprot_loading)

## joining gene and prot

library(dplyr)
gene_cprot = inner_join(gene2,cprot2, by = "name",multiple = "all")

##########
#finding pearson's correlation

# initialize two empty vectors to store correlation coefficient and p-value
corr <- vector()
pval <- vector()


sum(df1$gene_loading,na.rm = TRUE)
for (i in 1:nrow(df1)){
  test <- cor.test(df1$gene_loading[i], df1$cprot_loading[i])
  #extract correlation coefficient and p-value from test
  corr[i] = test$estimate
  pval[i] = test$p.value
}

?cor.test()
for (i in 1:nrow(gene_cprot)){
  test <- cor.test(gene_cprot$gene_loading[i], gene_cprot$cprot_loading[i])
  #extract correlation coefficient and p-value from test
  corr[i] = test$estimate
  pval[i] = test$p.value
}

#create new df with gene_cprot, corr, p_val

gene_cprot_corr = cbind(gene_cprot,corr, pval)







##########
#finding pearson's correlation

# initialize two empty vectors to store correlation coefficient and p-value
corr <- vector()
pval <- vector()

# loop over each element of b and c
for (i in 1:length(gene_cprot$gene_loading)) {
  
  # apply cor.test() function to b[i] and c[i] after converting them to numeric vectors
  test <- cor.test(as.numeric(unlist(gene_cprot$gene_loading[i])), as.numeric(unlist(gene_cprot$cprot_loading[i])))
  
  
  #extract correlation coefficient and p-value from test
  corr[i] = test$estimate
  pval[i] = test$p.value
}

#create new df with gene_cprot, corr, p_val

gene_cprot_corr = cbind(gene_cprot,corr, pval)

sum(gene_cprot, na.rm = TRUE)
