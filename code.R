#loading in all the libraries
library(BiocManager)

#libraries for DEG analysis
BiocManager::install("DESeq2")
BiocManager::install("GenomicFeatures")
BiocManager::install("tximport")
BiocManager::install("ComplexHeatmap")

#heatmap
library(ComplexHeatmap)
library(circlize)

#plots
install.packages("ggrepel")
library(ggrepel)

library(RColorBrewer)
library(tximport)
library(DESeq2)
library(GenomicFeatures) #annotating genes
library(ggplot2) #for pca and boxplots


install.packages("ggpubr")
library(ggpubr)

#loading in counts data
counts <- read.table("counts_I56.txt", sep = "\t", header = TRUE, row.names = 1)
dim(counts)

#editing counts data to remove irrelevant lines
counts2 <- counts[-c(21937, 21936, 21935),]
dim(counts2)

#loading in metadata
metadata <- read.table("metadata.csv", sep = ",", header = TRUE)

#setting the target and reference
reference <- "KGN"
target <- "COV434"

##PRE-FILTERING AND REMOVING GENES THAT ARE BARELY DETECTED
logic <- rowSums(counts2)> 10

#inspecting logical vector
summary(logic)
logic[1:3]

#subsetting the genes that are expressed
counts2 = counts2[logic,]

##NORMALIZING DATA USING DSEQ2

#creating DESEqDataset object
dds <- DESeqDataSetFromMatrix(countData = counts2, colData = metadata, design = ~ cell_line)

#examining dds
dds
dim(dds)
class(dds)
colData(dds)

#normalization and DESeq2 model fitting
dds = DESeq(dds)

res <- results (dds, contrast = c ("cell_line", target, reference))
resultsNames(dds)
res

#ordering by smallest pval
resOrdered <- res[order(res$pvalue),]
resOrdered
resOrdered[1:1,]

#checking number of adj pval lower than 0.05 threshold
sum(res$padj < 0.05, na.rm = TRUE)

#setting the alpha to 0.5
res05 <- results(dds, alpha = 0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm = TRUE)

res05 <- res[res$padj < 0.05,]

write.csv(res05, file = "deseq_results.csv") #saving the deseq results

##PCA ANALYSIS FOR QUALITY CONTROL

#gene expression matrix
vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = c("sample", "cell_line"))
pca_plot 

#PCA by prcomp
exprs.data = assay(vsd)
exprs.data[1:3,1:6]
pca_prcomp <- prcomp(t(exprs.data), scale = FALSE)

#inspecting prcomp
summary(pca_prcomp$x)
summary(pca_prcomp$rotation)
summary(pca_prcomp$sdev)
pca_prcomp$sdev[1:6]

sample_colors <- rainbow(length(unique(colData(vsd)$sample)))

variance_pc1 <- round(summary(pca_prcomp)$importance[2, 1], 2)
variance_pc2 <- round(summary(pca_prcomp)$importance[2, 2], 2)

pdf("PCA.pdf")
par(mfrow = c(2,2))
plot(pca_prcomp$x[,1:2], main = 'Sample Plot', col = sample_colors[as.factor(colData(vsd)$sample)], pch = 19)
text(pca_prcomp$x[,1:2], labels = rownames(pca_prcomp$x[,1:2]), pos = 4, cex = 0.8)
abline (v = 0, h = 0, col = 8)


plot(pca_prcomp$rotation [ ,1:2], main = " Gene Loading ")
abline ( v =0 , h =0 , col =8)

eig = (pca_prcomp$sdev)^2; names(eig) = 1: length(eig); eig = 100* eig/sum (eig)
barplot(eig, main = "Scree Plot", xlab = "PC", ylab = "Percentage of explained variances", ylim = c(0, max(eig) + 10))

dev.off()

##heat map
df <- as.data.frame(res05)
df

df.top <- df[(df$baseMean > 50) & (abs(df$log2FoldChange) > 2),]
df.top

dim(df.top)

df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]
df.top

rlog_out <- rlog(dds, blind = FALSE) #getting normalized count data from dds object
rlog_out

assay(rlog_out)

coldata <- data.frame(colnames(counts2))
coldata

mat <- assay(rlog_out)[rownames(df.top), coldata$colnames.counts2.] #sig genes x samples
mat
colnames(mat) <- coldata$colnames.counts2.
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
mat.scaled
colnames(mat.scaled)<-colnames(mat)

#keeping the top 25 most expressed genes and top 25 least expressed genes
num_keep<- 50
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)))
rows_keep

#the heatmap
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)

h<-h1
h


png("heatmap_v1.png", res = 300, width = 3000, height = 5500)
print(h)
dev.off()

##VOLCANO PLOT
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

pdf("keygenes_volcano_plot.pdf", width = 16, height = 12)

EnhancedVolcano(res05,
                lab = rownames(res05),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('AMH','CALB2','CHGA','FOXL2','INHBA','MUC1','NR5A1',
                              'SMARCA2','SMARCA4','SYP'),
                boxedLabels = TRUE,
                pointSize = 1.0,
                labSize = 6.0,
                col = c('black', 'black', 'black', 'red3'),
                colAlpha = 1) + coord_flip()
dev.off()

##INVESTIGATING GENES ONLY PRESENT IN 1 CELL LINE
normCounts <- counts(dds, normalized=TRUE)

#add kgn mean
kgnCol <- c("KGN_1", "KGN_2", "KGN_3")
kgnMean <- rowMeans(normCounts[, kgnCol, drop = FALSE])
normCounts <- cbind(normCounts, kgnMean)

#add cov434 mean
cov434Col <- c("COV434_1", "COV434_2", "COV434_3")
cov434Mean <- rowMeans(normCounts[, cov434Col, drop = FALSE])
normCounts <- cbind(normCounts, cov434Mean)

#filtering
sig_genes <- rownames(res05) #res05 has the DEGs with the significant padj
normCounts2 <- normCounts[sig_genes, , drop = FALSE]


lfc_genes <- df[(abs(df$log2FoldChange) > 2),]
sig_genes2 <- rownames(lfc_genes)
normCounts3 <- normCounts2[sig_genes2, , drop = FALSE] #this contains all the genes that have a lfc bigger than 2

filterKGN <- normCounts3[, "cov434Mean"] <= 0 #this filters out all the genes that have a mean of 0
KGNonly <- normCounts3[filterKGN, , drop = FALSE]
write.csv(KGNonly, file = "KGNonlyGenes.csv")

filterCOV434 <- normCounts3[, "kgnMean"] <= 0
COV434only <- normCounts3[filterCOV434, , drop = FALSE]
write.csv(COV434only, file = "COV434onlyGenes.csv")

##GSEA analysis
BiocManager::install("fgsea")
library(fgsea)
library(data.table)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

#generating list of genes with the deseq results exclusive to each cell line
kgnDeseqFilter <- rownames(KGNonly)
kgnDESEQres <- res05[kgnDeseqFilter, , drop = FALSE]
dim(kgnDESEQres)

write.csv(kgnDESEQres, file = "kgnDESEQres.csv")

cov434DeseqFilter <- rownames(COV434only)
cov434DESEQres <- res05[cov434DeseqFilter, , drop = FALSE]
dim(cov434DESEQres)

write.csv(cov434DESEQres, file = "cov434DESEQres.csv")

fgseaDeseq <- rbind(kgnDESEQres, cov434DESEQres)

#ranking genes by LFC
ranks <- fgseaDeseq$log2FoldChange
names(ranks) <- rownames(fgseaDeseq)

#performing hallmark GSEA
pathways <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt.txt")
fgseaRes <- fgseaMultilevel(pathways,  stats = ranks)
order <- order(fgseaRes$padj)
fgseaRes <- fgseaRes[order,]

#examining result
fgseaRes[1:3,]
fgseaRes$pathway <- sub(fgseaRes$pathway, pattern= 'HALLMARK_', replacement='')

fgseaRes$leadingEdge[8]

#visualizing pathways
fgsea_plot <- fgseaRes[abs(fgseaRes$NES) > 0.5,]
ord <- order(fgsea_plot$NES)
fgsea_plot <- fgsea_plot[ord,]

fp <- ggplot(fgsea_plot, mapping = aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(mapping = aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x = 'Pathway', y = 'Normalized Enrichment Score') +
  theme_minimal()

ggsave(filename = paste0("hallmark_gsea_plot.pdf"), plot = fp, device = "pdf", dpi = 300)

##MAKING MORE BOXPLOTS
muc1 <- "MUC1"

tmp <- plotCounts(dds, gene = muc1, intgroup = "cell_line", returnData = TRUE)
tmp$gene <- muc1

gene_plots <- ggplot(data = tmp, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(muc1)

ggsave(filename = paste0("gene_expr_boxplots_", muc1, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)

chga <- "CHGA"

tmp <- plotCounts(dds, gene = chga, intgroup = "cell_line", returnData = TRUE)
tmp$gene <- chga

gene_plots <- ggplot(data = tmp, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(chga)

ggsave(filename = paste0("gene_expr_boxplots_", chga, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)


syp <- "SYP"

tmp <- plotCounts(dds, gene = syp, intgroup = "cell_line", returnData = TRUE)
tmp$gene <- syp

gene_plots <- ggplot(data = tmp, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(syp)

ggsave(filename = paste0("gene_expr_boxplots_", syp, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)


inhba <- "INHBA"

tmp <- plotCounts(dds, gene = inhba, intgroup = "cell_line", returnData = TRUE)
tmp$gene <- inha

gene_plots <- ggplot(data = tmp, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(inhba)

ggsave(filename = paste0("gene_expr_boxplots_", inhba, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)

calb2 <- "CALB2"

tmp <- plotCounts(dds, gene = calb2, intgroup = "cell_line", returnData = TRUE)
tmp$gene <- calb2

gene_plots <- ggplot(data = tmp, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(calb2)

ggsave(filename = paste0("gene_expr_boxplots_", calb2, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)

sf1 <- "NR5A1"

tmp <- plotCounts(dds, gene = sf1, intgroup = "cell_line", returnData = TRUE)
tmp$gene <- sf1

gene_plots <- ggplot(data = tmp, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(sf1)

ggsave(filename = paste0("gene_expr_boxplots_", sf1, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)


amh <- "AMH"

tmp <- plotCounts(dds, gene = amh, intgroup = "cell_line", returnData = TRUE)
tmp$gene <- amh

gene_plots <- ggplot(data = tmp, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(amh)

ggsave(filename = paste0("gene_expr_boxplots_", amh, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)

foxl2 <- "FOXL2"

tmp <- plotCounts(dds, gene = foxl2, intgroup = "cell_line", returnData = TRUE)
tmp$gene <- foxl2

gene_plots <- ggplot(data = tmp, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(foxl2)

ggsave(filename = paste0("gene_expr_boxplots_", foxl2, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)

smarca2 <- plotCounts(dds, gene = smarca2, intgroup = "cell_line", returnData = TRUE)
smarca2$gene <- smarca2

gene_plots <- ggplot(data = smarca2, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(smarca2)

ggsave(filename = paste0("gene_expr_boxplots_", smarca2, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)

smarca4 <- plotCounts(dds, gene = smarca4, intgroup = "cell_line", returnData = TRUE)
smarca4$gene <- smarca4

gene_plots <- ggplot(data = smarca4, aes(x = cell_line, y = count, fill = cell_line)) +
  geom_boxplot(
    size = 0.25,
    outlier.size = 0.75,
    outlier.alpha = 0.75
  ) +
  scale_y_log10() +
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "",
    y = "Normalised count"
  ) +
  theme_bw() +
  ggtitle(smarca4)

ggsave(filename = paste0("gene_expr_boxplots_", smarca4, ".pdf"), plot = gene_plots, device = "pdf", dpi = 300)
