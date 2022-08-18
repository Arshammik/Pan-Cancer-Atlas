cnt <- read.delim("~/desktop/Systematic Review/TCGA-LUAD/Results/expression(log2+1)(cnt).csv")
c <- read.delim("~/desktop/Systematic Review/TCGA-LUAD/Results/c (tumor type and ID).csv")
sorted_dif <- read.delim("~/desktop/Systematic Review/TCGA-LUAD/Results/sorted_dif.csv")



setwd("~/desktop/Systematic Review/TCGA-LUAD/Results/")


pdf("boxplotwide.pdf", width = 100, height = 25)
boxplot(cnt)
dev.off()

# PCA for genes
pc <- prcomp(cnt)
ex.scale <- t(scale(t(cnt), scale = F))
pc1  <- prcomp(ex.scale)

pdf("PC and ex.scale.pdf", width = 15, height = 15)
plot(pc)
plot(ex.scale)
plot(pc1)
dev.off()

# PCA for samples
gr <- c$Sample.Type
grt <- paste0(c$Sample.Type, ".", seq_along(c$Sample.Type))
library(sure)
PCR  <- data.frame(pc1$r[,1:3], group = gr)
PCRT <- data.frame(pc1$r[,1:3], group = grt)
pdf("0.merged.pdf", width = 100, height = 100)
p1 <- ggplot(PCR, aes(PC1,PC2, color = group)) + geom_point(size = 9)  
p2 <- ggplot(PCRT, aes(PC1,PC2, label = group)) + geom_text(size = 3) 
grid.arrange(p1,p2, nrow = 2)
dev.off()

#Valcano Plot

sorted_dif$express <- "NO"
sorted_dif$express [sorted_dif$log2FoldChange > 1  & sorted_dif$padj < 0.05] <- "UP"
sorted_dif$express [sorted_dif$log2FoldChange < -1 & sorted_dif$padj < 0.05] <- "DOWN"

pdf("valcanoplot normal.pdf", width = 15, height = 15)
ggplot(sorted_dif, aes(log2FoldChange, -log10(padj), col=express)) +
  geom_point() + 
  theme_minimal()+ geom_vline(xintercept = c(1,-1), col = "red")+
  geom_hline(yintercept = -log10(0.05), col = "red") + guides(colour = guide_legend(override.aes = list(size=1.5))) 

dev.off()

pdf("valcanoplot colorful.pdf", width = 15, height = 15)
diseased_vs_healthy <- sorted_dif %>%
  mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                               log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                               TRUE ~ "ns")) 


ggplot(sorted_dif, aes(log2FoldChange, -log10(padj))) +
  geom_point() + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)))

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

diseased_vs_healthy %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                     limits = c(-10, 10))  
dev.off()



