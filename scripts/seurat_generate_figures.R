# 为 07-seurat.html 生成示例图，输出到 docs/images/seurat/
# 在项目根目录运行：Rscript scripts/seurat_generate_figures.R

options(warn = 1)
proj_root <- if (dir.exists("docs")) getwd() else "D:/GoogleDrive/GithubDesktopRepo/Bioinformatics"
out_dir <- file.path(proj_root, "docs", "images", "seurat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", out_dir)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
})

# 若未安装 pbmc3k 则先下载；从 pbmc3k.SeuratData 包加载
if (!requireNamespace("pbmc3k.SeuratData", quietly = TRUE)) {
  message("Downloading pbmc3k dataset...")
  InstallData("pbmc3k")
}
data("pbmc3k", package = "pbmc3k.SeuratData")
pbmc <- pbmc3k
# 转为当前 Seurat 版本格式（若为旧对象）
tryCatch({
  pbmc <- UpdateSeuratObject(pbmc)
}, error = function(e) NULL)
if (!"percent.mt" %in% colnames(pbmc[[]])) {
  tryCatch({
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  }, error = function(e) {
    message("Skipping percent.mt (using 0): ", conditionMessage(e))
    pbmc[["percent.mt"]] <- 0
  })
}

# 6.2.1 质控小提琴图（用 ggplot2 直接画 meta，避免 Seurat VlnPlot 与 ggplot2 兼容性问题）
md <- as.data.frame(pbmc[[]])
p1a <- ggplot(md, aes(x = 1, y = .data[["nFeature_RNA"]])) + geom_violin(fill = "steelblue") + theme_minimal() + labs(title = "nFeature_RNA") + xlab("")
p1b <- ggplot(md, aes(x = 1, y = .data[["nCount_RNA"]])) + geom_violin(fill = "steelblue") + theme_minimal() + labs(title = "nCount_RNA") + xlab("")
p1c <- ggplot(md, aes(x = 1, y = .data[["percent.mt"]])) + geom_violin(fill = "steelblue") + theme_minimal() + labs(title = "percent.mt") + xlab("")
if (requireNamespace("cowplot", quietly = TRUE)) {
  p1 <- cowplot::plot_grid(p1a, p1b, p1c, ncol = 3)
} else {
  p1 <- p1a
}
ggsave(file.path(out_dir, "01-qc-vln.png"), plot = p1, width = 10, height = 4, dpi = 150)
message("Saved 01-qc-vln.png")

# 6.2.2 质控散点图
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(file.path(out_dir, "02-qc-scatter.png"), plot = p2, width = 5, height = 4.5, dpi = 150)
message("Saved 02-qc-scatter.png")

# 标准流程
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# 6.4.1 分群图
p3 <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters")
ggsave(file.path(out_dir, "03-dim-umap.png"), plot = p3, width = 6, height = 5, dpi = 150)
message("Saved 03-dim-umap.png")

# 6.4.2 FeaturePlot
p4 <- FeaturePlot(pbmc, features = c("CD3D", "CD14", "MS4A1"))
ggsave(file.path(out_dir, "04-feature-plot.png"), plot = p4, width = 12, height = 4, dpi = 150)
message("Saved 04-feature-plot.png")

# 6.4.3 VlnPlot CD3D（用 ggplot2 避免 Seurat VlnPlot 兼容性）
cd3d <- as.numeric(GetAssayData(pbmc, slot = "data")["CD3D", ])
md2 <- cbind(pbmc[[]], CD3D = cd3d)
p5 <- ggplot(md2, aes(x = factor(seurat_clusters), y = CD3D, fill = factor(seurat_clusters))) + geom_violin(scale = "width") + theme_minimal() + labs(x = "Cluster", title = "CD3D") + theme(legend.position = "none")
ggsave(file.path(out_dir, "05-vln-cd3d.png"), plot = p5, width = 6, height = 4, dpi = 150)
message("Saved 05-vln-cd3d.png")

# 6.4.4 RidgePlot（用 ggplot2 脊线）
if (requireNamespace("ggridges", quietly = TRUE)) {
  p6 <- ggplot(md2, aes(x = CD3D, y = factor(seurat_clusters), fill = factor(seurat_clusters))) + ggridges::geom_density_ridges() + theme_minimal() + labs(x = "CD3D", y = "Cluster", title = "CD3D") + theme(legend.position = "none")
} else {
  p6 <- p5
}
ggsave(file.path(out_dir, "06-ridge-cd3d.png"), plot = p6, width = 6, height = 5, dpi = 150)
message("Saved 06-ridge-cd3d.png")

# 6.5 Marker
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers[order(markers$cluster, -markers$avg_log2FC), ]
top10 <- do.call(rbind, lapply(split(top10, top10$cluster), head, 10))
top_genes <- unique(top10$gene)

# 6.5.3 DotPlot
p7 <- DotPlot(pbmc, features = top_genes) + RotatedAxis()
ggsave(file.path(out_dir, "07-dotplot.png"), plot = p7, width = 10, height = 6, dpi = 150)
message("Saved 07-dotplot.png")

# 6.5.4 DoHeatmap
p8 <- DoHeatmap(pbmc, features = top_genes, size = 3)
ggsave(file.path(out_dir, "08-heatmap.png"), plot = p8, width = 10, height = 8, dpi = 150)
message("Saved 08-heatmap.png")

message("All figures saved to ", out_dir)
