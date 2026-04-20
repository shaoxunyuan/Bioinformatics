counts <- read.table("sample_rnaseq_counts.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(counts) <- counts$Gene_ID
counts$Gene_ID <- NULL

coldata <- read.table("sample_rnaseq_metadata.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
coldata$Condition <- factor(coldata$Condition, levels = c("Control", "Treatment"))

stopifnot(all(colnames(counts) == rownames(coldata)))

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)), colData = coldata, design = ~ Condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Condition", "Treatment", "Control"))
res <- as.data.frame(res)
res$Gene_ID <- rownames(res)
res <- res[order(res$padj, na.last = TRUE), ]

write.table(res, file = "03_deseq2_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

res$neglog10padj <- -log10(res$padj)
res$neglog10padj[!is.finite(res$neglog10padj)] <- NA

p <- ggplot(res, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(alpha = 0.7, size = 2, na.rm = TRUE) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Volcano plot (DESeq2)",
    x = "log2 Fold Change (Treatment vs Control)",
    y = "-log10(adjusted p-value)"
  )

ggsave("03_volcano.png", p, width = 7, height = 5, dpi = 200)
