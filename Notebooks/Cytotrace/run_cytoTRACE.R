df08 <- read.table("st08counts_cyto.txt", sep = "\t",header =T )
rownames(df08) <- read.table("st08genes_cyto.txt", sep = "\t",header =F )$V1

df10 <- read.table("st10.5counts_cyto.txt", sep = "\t",header =T )
rownames(df10) <- read.table("st10.5genes_cyto.txt", sep = "\t",header =F )$V1

df12 <- read.table("st12counts_cyto.txt", sep = "\t",header =T )
rownames(df12) <- read.table("st12genes_cyto.txt", sep = "\t",header =F )$V1

df13 <- read.table("st13counts_cyto.txt", sep = "\t",header =T )
rownames(df13) <- read.table("st13genes_cyto.txt", sep = "\t",header =F )$V1

df16 <- read.table("st16counts_cyto.txt", sep = "\t",header =T )
rownames(df16) <- read.table("st16genes_cyto.txt", sep = "\t",header =F )$V1

df18 <- read.table("st18counts_cyto.txt", sep = "\t",header =T )
rownames(df18) <- read.table("st18genes_cyto.txt", sep = "\t",header =F )$V1

df20 <- read.table("st20counts_cyto.txt", sep = "\t",header =T )
rownames(df20) <- read.table("st20genes_cyto.txt", sep = "\t",header =F )$V1

df22 <- read.table("st22counts_cyto.txt", sep = "\t",header =T )
rownames(df22) <- read.table("st22genes_cyto.txt", sep = "\t",header =F )$V1

df24 <- read.table("st24counts_cyto.txt", sep = "\t",header =T )
rownames(df24) <- read.table("st24genes_cyto.txt", sep = "\t",header =F )$V1

df27 <- read.table("st27counts_cyto.txt", sep = "\t",header =T )
rownames(df27) <- read.table("st27genes_cyto.txt", sep = "\t",header =F )$V1

Sys.setenv(RETICULATE_PYTHON="/anaconda3/bin/python")
library(CytoTRACE)



datasets <- list(df08, df10, df12, df13, df16, df18, df20, df22, df24, df27)
results <- iCytoTRACE(datasets)

plotCytoTRACE(results)
