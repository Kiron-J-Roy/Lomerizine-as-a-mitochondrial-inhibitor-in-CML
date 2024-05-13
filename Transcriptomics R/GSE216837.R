#Packages
library(DESeq2)
library(limma) 
library(tidyverse)
library(stringr)
library(ggplot2)
library(dplyr)
library(apeglm)

#read in gene count matrix
countTable <- read.csv('raw gene count matrix.csv')

#removing chr numbers from the end of ensemble gene IDs and making gene IDs row names 
countTable$gene_id <- str_pad(countTable$gene_id,20, "right", pad = "_")
countTable$gene_id <- str_sub(countTable$gene_id, end=-6)
countTable <- as.data.frame(countTable)
row.names(countTable) <- countTable$gene_id
countTable[,1] <- NULL

#creating experimental design table
id <- colnames(countTable)
gr <- c('Untreated','Treated','Untreated','Treated',
        'Untreated','Treated','Untreated','Treated')

colTable <-data.frame('ID' = id, 'Group' = gr,row.names=1)

##Making deseq object
dds_raw <- DESeqDataSetFromMatrix(countData = countTable, colData = colTable, design= ~ Group)

#Filtering 0-count genes
notAllZero <- rowSums(counts(dds_raw))>0
dds_raw <- dds_raw[notAllZero,]
dds_raw <- DESeq(dds_raw)

#Dispersion plot
plotDispEsts(dds_raw, main="Gene count dispersion")

#Plotting count density
rld_raw <-rlog(dds_raw)
plotDensities(assay(rld_raw),group=rld_raw$Group,legend='topright',main="Density of normalised gene counts")

#Correlation heatmap
cormat1 <- round(cor(countTable[,c(1:8)]),2)
melted_cormat1 <- melt(cormat1)
cor1 <- ggplot(data = melted_cormat1, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +labs(x= "Samples", y="Samples")
cor1 + theme(axis.text.x=element_text(angle=90, hjust=1))

#PCA
p1 <- plotPCA(rld_raw,intgroup=c("Group"), returnData=TRUE) 
percentVar <- round(100 * attr(p1, "percentVar"))
p1 <- cbind(p1, pGroup)
names(p1)[6] <- "Patient number"

ggplot(p1, aes(PC1, PC2, colour=group)) + geom_point(aes(shape=`Patient number`), size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("Gene count PCA") +
  coord_fixed() + scale_colour_discrete(name="Treated/untreated",labels=c("Treated","Untreated")) 

##Swap levels to get the opposite direction of differential testing
res_shrink <- lfcShrink(dds_raw, coef="Group_Untreated_vs_Treated", type="apeglm")

##Downstream analysis