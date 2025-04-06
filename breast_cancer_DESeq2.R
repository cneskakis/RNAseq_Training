# Breast Cancer RNAseq DESeq2 Analysis

# Φόρτωση πακέτων
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("DESeq2", ask = FALSE)
install.packages("ggplot2")
install.packages("readr")

library(DESeq2)
library(ggplot2)
library(readr)

# Φόρτωση του counts αρχείου
counts <- read_csv("breast_cancer_counts_matrix.csv")
counts <- as.data.frame(counts)
rownames(counts) <- counts$X1
counts$X1 <- NULL

# Δήλωση δειγμάτων (5 Tumor, 5 Control)
condition <- factor(c(rep("Tumor", 5), rep("Control", 5)))
colData <- data.frame(row.names = colnames(counts), condition)

# DESeq2 αντικείμενο
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

# Φιλτράρισμα για χαμηλή έκφραση
dds <- dds[rowSums(counts(dds)) > 10, ]

# Τρέξιμο ανάλυσης
dds <- DESeq(dds)
res <- results(dds)

# Volcano plot
res$gene <- rownames(res)
res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Yes", "No")

ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot – Breast Cancer") +
  theme_minimal()

# PCA plot
vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)))
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], condition)

ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Samples")