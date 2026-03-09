#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

message("== 01-alignment-demo.R ==")

# 定位项目根目录（默认：脚本位于 docs/，根目录为上一级）
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_path <- if (length(file_arg) == 1) normalizePath(file_arg, winslash = "/", mustWork = TRUE) else ""
script_dir <- if (nzchar(script_path)) dirname(script_path) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

if (!dir.exists(repo_root)) {
  stop("无法定位项目根目录：", repo_root, call. = FALSE)
}

out_dir <- file.path(repo_root, "docs", "images", "alignment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("输出目录：", out_dir)

pkgs_cran <- c("ape", "phangorn", "ggplot2")
pkgs_bioc <- c("Biostrings", "DECIPHER", "msa")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(msa)
  library(ape)
  library(phangorn)
  library(ggplot2)
})

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

set.seed(1)

## 1) 两两序列比对：DNA（全局/局部）
dna1 <- DNAString("ACGTTGACCTGAACTGAC")
dna2 <- DNAString("ACGTCGACCTGAACTGAC")

aln_global <- pairwiseAlignment(
  pattern = dna1,
  subject = dna2,
  type = "global",
  substitutionMatrix = nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE),
  gapOpening = -5,
  gapExtension = -2
)

aln_local <- pairwiseAlignment(
  pattern = dna1,
  subject = dna2,
  type = "local",
  substitutionMatrix = nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE),
  gapOpening = -5,
  gapExtension = -2
)

writeLines(c(
  "## Global alignment (DNA)",
  paste0("score = ", score(aln_global)),
  capture.output(aln_global),
  "",
  "## Local alignment (DNA)",
  paste0("score = ", score(aln_local)),
  capture.output(aln_local)
), file.path(out_dir, "pairwise_dna.txt"))

## 2) 两两序列比对：蛋白（BLOSUM62，局部）
prot1 <- AAString("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGE")
prot2 <- AAString("MKWVTFISLLFLFSSAYSRAIYEKLGLDYYKDDDDKSEVAH")
data("BLOSUM62", package = "Biostrings")

aln_prot <- pairwiseAlignment(
  pattern = prot1,
  subject = prot2,
  type = "local",
  substitutionMatrix = BLOSUM62,
  gapOpening = -10,
  gapExtension = -1
)

writeLines(c(
  "## Local alignment (Protein, BLOSUM62)",
  paste0("score = ", score(aln_prot)),
  capture.output(aln_prot)
), file.path(out_dir, "pairwise_protein.txt"))

## 3) 多序列比对：DNA（DECIPHER）
dna_set <- DNAStringSet(c(
  s1 = "ACGTTGACCTGAACTGAC",
  s2 = "ACGTCGACCTGAACTGAC",
  s3 = "ACGTTGACCTGAGCTGAC",
  s4 = "ACGTTGACCTGAACTAAC"
))

msa_dna <- AlignSeqs(dna_set, processors = 1)
aligned_dna <- msa_dna
writeXStringSet(aligned_dna, file.path(out_dir, "msa_dna.fasta"), format = "fasta")

## 4) 多序列比对：蛋白（msa + ClustalOmega）
aa_set <- AAStringSet(c(
  p1 = "MKWVTFISLLFLFSSAYSRA",
  p2 = "MKWVTFISLLFLFSSAYSRA",
  p3 = "MKWVTFISLLFLFSSAYARA",
  p4 = "MKWVTFISLLFLFSSAYSRV"
))

msa_aa <- msaClustalOmega(aa_set)
writeLines(capture.output(msa_aa), file.path(out_dir, "msa_protein.txt"))

## 5) 一致性（PID）矩阵
pid_matrix_from_aligned <- function(aln_strings, ignore_gaps = TRUE) {
  x <- do.call(rbind, strsplit(aln_strings, ""))
  n <- nrow(x)
  out <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(names(aln_strings), names(aln_strings)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      a <- x[i, ]
      b <- x[j, ]
      keep <- rep(TRUE, length(a))
      if (ignore_gaps) keep <- (a != "-" & b != "-")
      denom <- sum(keep)
      out[i, j] <- if (denom == 0) NA_real_ else 100 * sum(a[keep] == b[keep]) / denom
    }
  }
  out
}

aln_strings <- toupper(as.character(aligned_dna))
pid_dna <- pid_matrix_from_aligned(aln_strings, ignore_gaps = TRUE)
write.csv(round(pid_dna, 4), file.path(out_dir, "pid_matrix.csv"), row.names = TRUE)

## 6) 距离矩阵（K80）+ NJ 树
aln_for_phy <- gsub("-", "n", tolower(aln_strings)) # gap 视作未知碱基 n（并统一为小写）
m_char <- do.call(rbind, strsplit(aln_for_phy, ""))
rownames(m_char) <- names(aln_for_phy)
dnabin <- as.DNAbin(m_char)

dist_k80 <- dist.dna(dnabin, model = "K80", pairwise.deletion = TRUE)
write.csv(as.matrix(dist_k80), file.path(out_dir, "dist_k80.csv"), row.names = TRUE)

tree_nj <- nj(dist_k80)
png(file.path(out_dir, "tree_nj.png"), width = 980, height = 680, res = 130)
plot(tree_nj, main = "NJ tree (K80 distance)")
add.scale.bar()
dev.off()

## 7) 距离热图（ggplot2）
dm <- as.matrix(dist_k80)
df_dm <- as.data.frame(as.table(dm))
colnames(df_dm) <- c("seq1", "seq2", "dist")

p <- ggplot(df_dm, aes(x = seq1, y = seq2, fill = dist)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient(low = "#0ea5e9", high = "#f97316") +
  coord_equal() +
  labs(title = "K80 distance matrix", x = NULL, y = NULL, fill = "distance") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(out_dir, "dist_heatmap.png"), p, width = 7.6, height = 6.2, dpi = 160)

## 8) 变异位点（SNP/Indel）表
m <- do.call(rbind, strsplit(aln_strings, ""))
col_is_var <- apply(m, 2, function(x) length(unique(x)) > 1)
col_has_gap <- apply(m, 2, function(x) any(x == "-"))
idx <- which(col_is_var)

var_tbl <- data.frame(
  pos = idx,
  has_gap = col_has_gap[idx],
  pattern = apply(m[, idx, drop = FALSE], 2, function(x) paste(x, collapse = "")),
  stringsAsFactors = FALSE
)
write.csv(var_tbl, file.path(out_dir, "variant_sites.csv"), row.names = FALSE)

## 9) Motif 检索（在原始序列中匹配）
motif <- DNAString("ACCTG")
motif_hits <- do.call(rbind, lapply(names(dna_set), function(nm) {
  x <- dna_set[[nm]]
  h <- matchPattern(motif, x)
  if (length(h) == 0) {
    return(data.frame(seq = nm, start = integer(), width = integer()))
  }
  data.frame(seq = nm, start = start(h), width = width(h))
}))
write.csv(motif_hits, file.path(out_dir, "motif_hits.csv"), row.names = FALSE)

message("完成：已生成示例输出文件：")
message(" - ", file.path(out_dir, "pairwise_dna.txt"))
message(" - ", file.path(out_dir, "pairwise_protein.txt"))
message(" - ", file.path(out_dir, "msa_dna.fasta"))
message(" - ", file.path(out_dir, "pid_matrix.csv"))
message(" - ", file.path(out_dir, "dist_k80.csv"))
message(" - ", file.path(out_dir, "dist_heatmap.png"))
message(" - ", file.path(out_dir, "tree_nj.png"))
message(" - ", file.path(out_dir, "variant_sites.csv"))
message(" - ", file.path(out_dir, "motif_hits.csv"))

