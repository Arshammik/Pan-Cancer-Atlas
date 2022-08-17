################################################################################
################################################################################
setRepositories()

library(data.table)
library(limma)
library(VennDiagram)
library(ggplot2)
library(parallel)
library(pheatmap)
library(reshape2)
library(dplyr)


#TUAD
#options(stringsAsFactors = F)
setwd("~/desktop/Systematic Review/")
ex <- read.delim("~/Desktop/Systematic Review/ex_LUAD.tsv")
counts <- ex


#counts3 <- aggregate(.~gene_name,counts,mean)
# Can't run on local

dim(counts)
counts <- distinct(counts, gene_name, .keep_all = TRUE)
dim(counts)
# lose 1233 genes

rownames(counts) <- counts$gene_name
counts <- counts[,-1:-2]

clinical <- read.delim("../Systematic Review/TCGA-LUAD/Data/clinical.cart.2022-06-11/clinical.tsv")
exposure <- read.delim("../Systematic Review/TCGA-LUAD/Data/clinical.cart.2022-06-11/exposure.tsv")
family_history <- read.delim("../Systematic Review/TCGA-LUAD/Data/clinical.cart.2022-06-11/family_history.tsv")
sampleID <- read.delim("../Systematic Review/TCGA-LUAD/Data/gdc_sample_sheet.2022-06-11.tsv")
samples <- sampleID[,7:8]

a <- colnames(counts)
a <- sub('\\.', '-', a)
a <- sub('\\.', '-', a)
a <- sub('\\.', '-', a)
a <- sub('\\.', '-', a)

b <- data.frame(Sample.ID = a)

samples <- distinct(samples, Sample.ID, .keep_all = TRUE)
c <- merge(b, samples, by.x = "Sample.ID", by.y = "Sample.ID", sort = F)


colnames (counts) <- c$Sample.Type
head(counts)

# Data Prepration is ready
## Next step will be performed by DESeq2 library