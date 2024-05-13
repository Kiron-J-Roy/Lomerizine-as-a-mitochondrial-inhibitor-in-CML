### 1 install mEndtoEnd packages ###
BiocManager::install("maEndToEnd")
liberary("maEndToEnd")

### 2 List of packages required for the workflow ###

#General Bioconductor packages
library(Biobase)
library(oligoClasses)
#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)
#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)
#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
#Formatting/documentation packages
library(rmarkdown)
library(BiocStyle)
library(dplyr)
library(tidyr)
#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)

### 3 Downloading the raw data from ArrayExpress ###

#setwd
raw_data_dir=setwd("C:/Users/user1/Desktop/HSC_CML_TKI/preprocessing of arrays/HSC_Array2")
getwd()

#get raw files
library(ArrayExpress)
anno_AE <- getAE("E-MTAB-2581", path = raw_data_dir, type = "raw")
#get SDRF that contains info about samples
sdrf_location <- file.path(raw_data_dir, "E-MTAB-2581.sdrf.txt")
SDRF <- read.delim(sdrf_location)
rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)

### 4 read cel files and generate expressionset as raw_data ### 
raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir,
                                                       SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)

#=============Warning message:In oligo::read.celfiles(filenames = file.path(raw_data_dir, SDRF$Array.Data.File),:'channel' automatically added to varMetadata in phenoData. it means we have to select columns we are interested in and edit/replace with shorted names as follows:

View(SDRF@varMetadata) # view which pdata you are interested in
View(raw_data@phenoData@data) # view experimental samples conditions
#filter to required columns 
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[ ,c("Factor.Value.cell.type.","Factor.Value.disease.")]
#rename phnodata of raw_data into pdata
pdata=raw_data@phenoData@data
view(pdata)
colnames(pdata)=c("cell","disease")
#writepdata table and edit within excel to simplistic version
write.table(pdata,file="pdata2.csv", quote = FALSE, sep = ",")
#======= do not forget to edit with excel before excute next =========#
#read it again
pdata2=read.table("pdata2.csv", head=TRUE, row.name=1, sep = ",")
#replace pdata with pdata2
pdata=pdata2 
raw_data@phenoData@data=pdata2
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[ ,c("cell","disease")]
################ Quality control of the raw data #########################
##########################################################################
# log2 QC then PCA 
exp_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     cell = pData(raw_data)$cell,
                     disease= pData(raw_data)$disease)
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = cell, colour = disease), size=5) +
  ggtitle("raw expression") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  scale_shape_manual(values = c(0,15)) +
  scale_color_manual(values = c("brown2", "dodgerblue4"))

#boxplot of arrays expression intensities
oligo::boxplot(raw_data, target = "core",
               main = "intensities raw data")

#### Report of raw_eset #####
arrayQualityMetrics(expressionset = raw_data,
                    outdir = setwd("C:/Users/user1/Desktop/HSC_Array2"),
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("cell", "disease"))

#Background correction without normalization of expression 
gr_eset <- oligo::rma(raw_data, target = "core", normalize = FALSE)

# boxplot relative expression deviation for Norm_eset
row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(gr_eset)))
RLE_data <- sweep(Biobase::exprs(gr_eset), 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- tidyr::gather(RLE_data, patient_array, log2_expression_deviation)
ggplot2::ggplot(RLE_data_gathered, aes(patient_array,log2_expression_deviation)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(c(-2, 2)) +
  theme(axis.text.x = element_text(colour = "aquamarine4",angle = 60, size = 6.5, hjust = 1 ,face = "bold"))

### 6. One-step preprocessing normalization in oligo ###
rma <- oligo::rma(raw_data, target = "core") #this will be used for limma after
rma2=rma@assayData[["exprs"]]
write.csv(rma2,"expression.csv")
# log2 QC then PCA 
exp_rma <- log2(Biobase::exprs(rma))
PCA_rma <- prcomp(t(exp_rma), scale. = FALSE)
percentVar2 <- round(100*PCA_rma$sdev^2/sum(PCA_rma$sdev^2),1)
sd_ratio2 <- sqrt(percentVar2[2] / percentVar2[1])
dataGG2 <- data.frame(PC1 = PCA_rma$x[,1], 
                      PC2 = PCA_rma$x[,2],
                     cell = pData(rma)$cell,
                     disease= pData(rma)$disease)
ggplot(dataGG2, aes(PC1, PC2)) +
  geom_point(aes(shape = cell, colour = disease), size=5) +
  ggtitle("rma expression") +
  xlab(paste0("PC1, VarExp: ", percentVar2[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar2[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  scale_shape_manual(values = c(0,15)) +
  scale_color_manual(values = c("brown2", "dodgerblue4"))


#boxplot of arrays expression intensities
oligo::boxplot(rma, target = "core",
               main = "intensities raw data")

#### Report of raw_eset #####
arrayQualityMetrics(expressionset = rma,
                    outdir = setwd("C:/Users/user1/Desktop/HSC_Array2"),
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("cell", "disease"))
getwd()


#RLE after normalization
row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(rma)))
RLE_data <- sweep(Biobase::exprs(rma), 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- tidyr::gather(RLE_data, patient_array, log2_expression_deviation)
library(ggplot2)
ggplot(RLE_data_gathered, aes(patient_array,log2_expression_deviation)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(c(-2, 2)) +
  theme(axis.text.x = element_text(colour = "aquamarine4",angle = 60, size = 6.5, hjust = 1 ,face = "bold"))




########## limma #################
# subset then linear modeling of SC
subset.rma <- rma[, rma$cell %in% c("SC") & rma$disease %in% c("CML","normal")]

#subset for SCvsLPC comparison
subset.rma2 <- rma[, rma$cell %in% c("SC","PC") & rma$disease %in% c("CML")]

# generate matrix for cell type (HSC vs LSC) (FC is HSC-LSC)
design <- model.matrix(~ subset.rma$disease)
head(design)
fit <- lmFit(subset.rma, design)
fit <- eBayes(fit)
DE=topTable(fit,sort="none",n=Inf)

write.table(DE,"HSCvsCML.csv")
HSCvsCML=read.table("HSCvsCML.csv")

#for that reason reverse matrix should be done (LSC vs HSC)
design2 <- model.matrix(~ subset.rma$disease-1)
head(design2)
colnames(design2)=c("LSC","HSC")
design2
fit2 <- lmFit(subset.rma, design2)
contrast.matrix <- makeContrasts("LSC-HSC", levels = design2)
contrast.matrix
fit2C <- contrasts.fit(fit2, contrast.matrix)
fit2C <- eBayes(fit2C)
DE2=topTable(fit2C,sort="none",n=Inf)

write.table(DE2,"LSCvsHSC.csv")
LSCvsHSC=read.table("LSCvsHSC.csv")
LSCvsHSC

# another everse matrix should be done (LSC vs LPC)
design3 <- model.matrix(~ subset.rma2$cell-1)
head(design3)
colnames(design3)=c("LSC","LPC")
design3
fit3 <- lmFit(subset.rma2, design3)
contrast.matrix2 <- makeContrasts("LSC-LPC", levels = design3)
contrast.matrix2
fit3C <- contrasts.fit(fit3, contrast.matrix2)
fit3C <- eBayes(fit3C)
DE3=topTable(fit3C,sort="none",n=Inf)

write.table(DE3,"LSCvsLPC.csv")
LSCvsLPC=read.table("LSCvsLPC.csv")
LSCvsLPC

# subset then linear modeling of PC with reverse matrix
subset.rma3 <- rma[, rma$cell %in% c("PC") & rma$disease %in% c("CML","normal")]

#for that reason reverse matrix should be done (LPC vs HPC)
design4 <- model.matrix(~ subset.rma4$disease-1)
head(design4)
colnames(design4)=c("LPC","HPC")
design4
fit4 <- lmFit(subset.rma4, design4)
contrast.matrix4 <- makeContrasts("LPC-HPC", levels = design4)
contrast.matrix4
fitC4 <- contrasts.fit(fit4, contrast.matrix4)
fitC4 <- eBayes(fitC4)
DE4=topTable(fitC4,sort="none",n=Inf)

write.table(DE4,"LPCvsHPC.csv")
LPCvsHPC=read.table("LPCvsHPC.csv")
LPCvsHPC
#generate plots for foldchange and limma
## blot these results

#SA plot for LSC vs HSC
plotSA(fit)
ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
plotSA(fit)
#MD plot
plotMD(fit)
abline(0,0,col="blue")
#Venn diagram
results <- decideTests(fit)
vennDiagram(results)

#SA plot for LSC vs LPC
plotSA(fit3)
ordinary.t2 <- fit3$coef / fit3$stdev.unscaled / fit3$sigma
plotSA(fit3)
#MD plot for LSCvsLPC
plotMD(fit3)
abline(0,0,col="blue")
#Venn diagram
results3 <- decideTests(fit3)
vennDiagram(results3)

## 7. Filtering lowly expressed genes based on intensity ###
#expression intensities of all genes 
medians <- rowMedians(Biobase::exprs(rma))
histogram <- hist(medians, 100, col = "cornsilk1", freq = FALSE,
                  main = "median intensities",
                  border = "antiquewhite4",
                  xlab = "Median intensities")
# make cutoff point
man_threshold <- 4
histogram <- hist(medians, 100, col = "cornsilk", freq = FALSE,
                  main = "expression pattern",
                  border = "antiquewhite4",
                  xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)
# more trimming below threshold of sample number
# view no of samples to set threshold
no_of_samples <-
  table(paste0(pData(rma)$Factor.Value.compound.))
no_of_samples
# exclude values that do not reach that level (did not affect final rma outcome)
samples_cutoff <- min(no_of_samples)
threshold <- apply(Biobase::exprs(rma), 1,function(x){sum(x > man_threshold) >= samples_cutoff})
table(threshold)
# listen to Egypt scholar Mariam first lecture and see where she is doing it


############################# 8. annotation ############################
########################################################################

ID <- featureNames(rma)
Symbol <- getSYMBOL(ID, "hugene10sttranscriptcluster.db")
Name <- as.character(lookUp(ID, "hugene10sttranscriptcluster.db",
                            "GENENAME"))
#make a temporary data frame with all the identifiers...
tmpframe <-data.frame(ID=ID, Symbol=Symbol,
                      Name=Name,stringsAsFactors=F)
tmpframe[tmpframe=="NA"] <- NA
#assign data frame to rma-results
fData(rma) <- tmpframe
View(rma@featureData@data) # view new feature data added to rma assay data
View(raw_data@featureData) # nothing on it. 
rma.annotated = cbind(pData(featureData(rma))[,"Symbol"],exprs(rma)) # table
colnames(rma.annotated)[1] ="Symbol" # rename v1 with symbol gene
#generate table
write.table(rma.annotated,file="better_annotation.csv", quote = FALSE, sep = ",")

############################# 8. annotation ############################
########################################################################


ID <- featureNames(rma)
Symbol <- getSYMBOL(ID, "hugene10sttranscriptcluster.db")
Name <- as.character(lookUp(ID, "hugene10sttranscriptcluster.db",
                            "GENENAME"))
#make a temporary data frame with all the identifiers...
tmpframe <-data.frame(ID=ID, Symbol=Symbol,
                      Name=Name,stringsAsFactors=F)
tmpframe[tmpframe=="NA"] <- NA
#assign data frame to rma-results
fData(rma) <- tmpframe
View(rma@featureData@data) # view new feature data added to rma assay data
View(raw_data@featureData) # nothing on it. 
rma.annotated = cbind(pData(featureData(rma))[,"Symbol"],exprs(rma)) # table
colnames(rma.annotated)[1] ="Symbol" # rename v1 with symbol gene
#generate table
write.table(rma.annotated,file="better_annotation.csv", quote = FALSE, sep = ",")
