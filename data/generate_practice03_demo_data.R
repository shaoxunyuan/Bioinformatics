# 生成练习 03 配套示例：表达矩阵、差异表达结果、样本 metadata
# 运行：在项目根目录执行  Rscript data/generate_practice03_demo_data.R

set.seed(42)

out_dir <- if (basename(getwd()) == "data") dirname(getwd()) else getwd()
data_dir <- file.path(out_dir, "data")

genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")
cancer_means <- c(100.5, 12.3, 2.3, 3.1, 7.0)
health_means <- c(25.1, 45.7, 30.4, 3.0, 7.3)

n_cancer <- 4L
n_health <- 4L
mat <- matrix(NA_real_, nrow = length(genes), ncol = n_cancer + n_health)
rownames(mat) <- genes
colnames(mat) <- paste0("Sample", seq_len(n_cancer + n_health))

for (i in seq_along(genes)) {
  sc <- max(cancer_means[i] * 0.04, 0.2)
  sh <- max(health_means[i] * 0.04, 0.2)
  noise_c <- rnorm(n_cancer, sd = sc)
  noise_c <- noise_c - mean(noise_c)
  noise_h <- rnorm(n_health, sd = sh)
  noise_h <- noise_h - mean(noise_h)
  mat[i, seq_len(n_cancer)] <- pmax(0.01, cancer_means[i] + noise_c)
  mat[i, n_cancer + seq_len(n_health)] <- pmax(0.01, health_means[i] + noise_h)
}

mat <- round(mat, 2)
expr <- data.frame(Gene_ID = genes, mat, check.names = FALSE)

de <- data.frame(
  Gene_ID = genes,
  mean_Cancer = c(100.5, 12.3, 2.3, 3.1, 7.0),
  mean_Health = c(25.1, 45.7, 30.4, 3.0, 7.3),
  log2FoldChange = c(2.01, -1.89, 1.62, 0.97, 1.03),
  padj = c(1.2e-6, 4.5e-5, 0.11, 0.25, 3.9e-5),
  significance = c("Significant", "Significant", "Non-Significant", "Non-Significant", "Significant"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta <- data.frame(
  Sample = colnames(mat),
  Condition = c(rep("Cancer", n_cancer), rep("Health", n_health)),
  stringsAsFactors = FALSE
)

write.csv(expr, file.path(data_dir, "practice03_expression_matrix.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(de, file.path(data_dir, "practice03_de_results.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.table(
  meta,
  file.path(data_dir, "practice03_metadata.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)

message("Written:\n - ", file.path(data_dir, "practice03_expression_matrix.csv"), "\n - ",
        file.path(data_dir, "practice03_de_results.csv"), "\n - ",
        file.path(data_dir, "practice03_metadata.tsv"))
