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
cn <- ex


#cn3 <- aggregate(.~gene_name,cn,mean)
# Can't run on local

dim(cn)
cn <- distinct(cn, gene_name, .keep_all = TRUE)
dim(cn)
# lose 1233 genes

rownames(cn) <- cn$gene_name
cn <- cn[,-1:-2]

#clinical <- read.delim("../Systematic Review/TCGA-LUAD/Data/clinical.cart.2022-06-11/clinical.tsv")
#exposure <- read.delim("../Systematic Review/TCGA-LUAD/Data/clinical.cart.2022-06-11/exposure.tsv")
#family_history <- read.delim("../Systematic Review/TCGA-LUAD/Data/clinical.cart.2022-06-11/family_history.tsv")
sampleID <- read.delim("../Systematic Review/TCGA-LUAD/Data/gdc_sample_sheet.2022-06-11.tsv")
i1 <- sampleID$Data.Type == "Gene Expression Quantification"
samples <- sampleID[i1, , drop = F]
samples <- samples[,7:8]

a <- colnames(cn)
a <- sub('\\.', '-', a)
a <- sub('\\.', '-', a)
a <- sub('\\.', '-', a)
a <- sub('\\.', '-', a)

b <- data.frame(Sample.ID = a)
c <- merge(b, samples, by.x = "Sample.ID", by.y = "Sample.ID", sort = F)


colnames (cn) <- c$Sample.Type
colnames(cn)

# Data Prepration is ready
## Next step will be performed by DESeq2 library

dim(cn)
#install.packages("DESeq2")
library(DESeq2)

gr <- factor(c$Sample.Type)
colSums(cn) #Through this we can demonstrates 
colData <- data.frame(group = gr, type = "paired-end")
cds <- DESeqDataSetFromMatrix(cn, colData, design = ~group)
cds <- DESeq(cds)
cnt <- log2(1+(counts(cds, normalize = T))) #getting normalized counts

write.table(cnt, "~/desktop/Systematic Review/TCGA-LUAD/Results/expression(log2+1)(cnt).csv",  quote = F, col.names = T, row.names = T, sep = "\t")
write.table(c, "~/desktop/Systematic Review/TCGA-LUAD/Results/c (tumor type and ID).csv",  quote = F, col.names = T, row.names = T, sep = "\t")


##
#DEGs
dif <- results(cds, c("group", "Solid Tissue Normal" , "Primary Tumor"))
write.table(dif, "~/desktop/Systematic Review/TCGA-LUAD/Results/dif.csv",  quote = F, col.names = T, row.names = T, sep = "\t")

#checkpoint
sorted_dif <- data.frame(results(cds, c("group", "Solid Tissue Normal" , "Primary Tumor")))
sorted_dif$padj <- p.adjust(sorted_dif$pvalue, method = "BH")
sorted_dif <- sorted_dif[order(sorted_dif$padj),]

X  <- subset(sorted_dif, log2FoldChange > 1  & padj < 0.05)
XPrime  <- subset(dif, log2FoldChange > 1  & padj < 0.05)
Y  <- subset(sorted_dif, log2FoldChange < -1  & padj < 0.05)
YPrime  <- subset(dif, log2FoldChange < -1  & padj < 0.05)

dim(X); dim(XPrime); dim(Y); dim(YPrime) 

write.table(dif, "dif.csv", quote = F, col.names = T, row.names = T, sep = "\t" )
write.table(sorted_dif, "sorted_dif.csv", quote = F, col.names = T, row.names = T, sep = "\t" )

write.table(subset(X), "Up regulated(sorted_dif).csv",  quote = F, col.names = T, row.names = T, sep = "\t")
write.table(subset(XPrime), "Up regulated(dif).csv",  quote = F, col.names = T, row.names = T, sep = "\t")
write.table(subset(Y), "Down regulated(sorted_dif).csv",  quote = F, col.names = T, row.names = T, sep = "\t")
write.table(subset(YPrime), "Up regulated(dif).csv",  quote = F, col.names = T, row.names = T, sep = "\t")









