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
raw_data_dir=setwd("C:/Users/user1/Desktop/CML_Array")
getwd()

#get raw files
library(ArrayExpress)
anno_AE <- getAE("E-MTAB-2594", path = raw_data_dir, type = "raw")
#get SDRF that contains info about samples
sdrf_location <- file.path(raw_data_dir, "E-MTAB-2594.sdrf.txt")
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
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[ ,c("Factor.Value..cell.type.","Factor.Value.compound.","Factor.Value.time.")]
#rename phnodata of raw_data into pdata
pdata.ima=raw_data@phenoData@data
pdata.ima
colnames(pdata.ima)=c("cell","compound","time")
#writepdata table and edit within excel to simplistic version
write.table(pdata.ima,file="pdata.ima.csv", quote = FALSE, sep = ",")
#======= do not forget to edit with excel before excute next =========#
#read it again
pdata.ima=read.table("pdata.ima.csv", head=TRUE, row.name=1, sep = ",")
#replace pdata with pdata2
pdata=pdata.ima
raw_data@phenoData@data=pdata.ima
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[ ,c("cell","compound","time")]

################ Quality control of the raw data #########################
##########################################################################

#QC requires log then PCA on that (Error in ggplot(dataGG, aes(PC1, PC2)) : object 'dataGG' not found) could not regenerate it (do not know why)

# boxplot of relative expression deviation for Norm_eset
rma <- oligo::rma(raw_data, scale=FALSE)
row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(rma)))
RLE_data <- sweep(Biobase::exprs(rma), 1, row_medians_assayData)
RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- tidyr::gather(RLE_data, patient_array, log2_expression_deviation)
ggplot2::ggplot(RLE_data_gathered, aes(patient_array,log2_expression_deviation)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(c(-2, 2)) +
  theme(axis.text.x = element_text(colour = "aquamarine4",angle = 60, size = 6.5, hjust = 1 ,face = "bold"))

### 6. One-step preprocessing normalization in oligo ###
rma <- oligo::rma(raw_data, target = "core") #this will be used for limma after changing phenodata with pdata we generated before

# boxplot of relative expression deviation for rma_eset
rma_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(rma)))
RLE_data <- sweep(Biobase::exprs(rma), 1, rma_medians_assayData)
RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- tidyr::gather(RLE_data, patient_array, log2_expression_deviation)
library(ggplot2)
ggplot2::ggplot(RLE_data_gathered, aes(patient_array,log2_expression_deviation)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(c(-2, 2)) +
  theme(axis.text.x = element_text(colour = "aquamarine4",angle = 60, size = 6.5, hjust = 1 ,face = "bold"))

#boxplot of arrays expression intensities
oligo::boxplot(raw_data, target = "core",
               main = "intensities rma data")

#boxplot of arrays expression intensities
oligo::boxplot(rma, target = "core",
               main = "intensities rma data")



# subset then data as stem cells (SC) and imatinib 7 days
subset.rma <- rma[, rma$cell %in% c("SC") & rma$compound %in% c("imatinib","none")&rma$time%in%c("7","0")]
#for that reason reverse matrix should be done (LSC vs HSC)
design <- model.matrix(~ subset.rma$compound-1)
head(design)
colnames(design)=c("Ima7","LSC")
design
fit <- lmFit(subset.rma, design)
contrast.matrix <- makeContrasts("Ima7-LSC", levels = design)
contrast.matrix
fitC <- contrasts.fit(fit, contrast.matrix)
fitC <- eBayes(fitC)
DE2=topTable(fitC,sort="none",n=Inf)

write.table(DE2,"Ima7vsLSC.csv")
Ima7vsLSC=read.table("Ima7vsLSC.csv")
Ima7vsLSC

# subset then data as stem cells (SC) and nilotinib 7 days
subset.rma3 <- rma[, rma$cell %in% c("SC") & rma$compound %in% c("nilotinib","none")&rma$time%in%c("7","0")]
#for that reason reverse matrix should be done (LSC vs HSC)
design3 <- model.matrix(~ subset.rma3$compound-1)
head(design3)
colnames(design3)=c("nil","LSC")
design3
fit3 <- lmFit(subset.rma3, design3)
contrast.matrix3 <- makeContrasts("nil-LSC", levels = design3)
contrast.matrix3
fitC3 <- contrasts.fit(fit3, contrast.matrix3)
fitC3 <- eBayes(fitC3)
DE23=topTable(fitC3,sort="none",n=Inf)

write.table(DE23,"nil7vsLSC.csv")
nil7vsLSC=read.table("nil7vsLSC.csv")
nil7vsLSC


# subset then data as stem cells (SC) and nilotinib 8 hours
subset.rma4 <- rma[, rma$cell %in% c("SC") & rma$compound %in% c("nilotinib","none")&rma$time%in%c("8","0")]
#for that reason reverse matrix should be done (LSC vs HSC)
design4 <- model.matrix(~ subset.rma4$compound-1)
head(design4)
colnames(design4)=c("nil","LSC")
design4
fit4 <- lmFit(subset.rma4, design4)
contrast.matrix4 <- makeContrasts("nil-LSC", levels = design4)
contrast.matrix4
fitC4 <- contrasts.fit(fit4, contrast.matrix4)
fitC4 <- eBayes(fitC4)
DE24=topTable(fitC4,sort="none",n=Inf)

write.table(DE24,"nil8vsLSC.csv")
nil8vsLSC=read.table("nil8vsLSC.csv")
nil8vsLSC

# subset then data as stem cells (SC) and dasatinib 8 hours
subset.rma5 <- rma[, rma$cell %in% c("SC") & rma$compound %in% c("dasatinib","none")&rma$time%in%c("8","0")]
#for that reason reverse matrix should be done (LSC vs HSC)
design5 <- model.matrix(~ subset.rma5$compound-1)
head(design5)
colnames(design5)=c("das","LSC")
design5
fit5 <- lmFit(subset.rma5, design5)
contrast.matrix5 <- makeContrasts("das-LSC", levels = design5)
contrast.matrix5
fitC5 <- contrasts.fit(fit5, contrast.matrix5)
fitC5 <- eBayes(fitC5)
DE25=topTable(fitC5,sort="none",n=Inf)

write.table(DE25,"das8vsLSC.csv")
das8vsLSC=read.table("das8vsLSC.csv")
das8vsLSC

# subset then data as stem cells (SC) and dasatinib 7 days
subset.rma6 <- rma[, rma$cell %in% c("SC") & rma$compound %in% c("dasatinib","none")&rma$time%in%c("7","0")]
#for that reason reverse matrix should be done (LSC vs HSC)
design6 <- model.matrix(~ subset.rma6$compound-1)
head(design6)
colnames(design6)=c("das","LSC")
design6
fit6 <- lmFit(subset.rma6, design6)
contrast.matrix6 <- makeContrasts("das-LSC", levels = design6)
contrast.matrix6
fitC6 <- contrasts.fit(fit6, contrast.matrix6)
fitC6 <- eBayes(fitC6)
DE26=topTable(fitC6,sort="none",n=Inf)

write.table(DE26,"das7vsLSC.csv")
das7vsLSC=read.table("das7vsLSC.csv")
das7vsLSC


#generate plots for foldchange and limma
## blot these results

#SA plot
plotSA(fit)
ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
plotSA(fit)
#MD plot
plotMD(fit)
abline(0,0,col="blue")
#Venn diagram
results <- decideTests(fit)
vennDiagram(results)


#============== we already generated all statistical figures for plotting ================= 


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