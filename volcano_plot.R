# Volcano Plot in R using dummy RNAseq data

# Load libraries
library(ggplot2)

# Simulated RNAseq results
set.seed(123)
data <- data.frame(
  gene = paste0("Gene", 1:1000),
  log2FoldChange = rnorm(1000, mean = 0, sd = 2),
  pvalue = runif(1000, 0, 1)
)
data$significant <- ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) > 1, "Yes", "No")

# Volcano plot
ggplot(data, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(p-value)") +
  theme_minimal()
