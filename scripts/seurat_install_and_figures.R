# 若缺少 SeuratData 则安装，然后生成图
proj_root <- "D:/GoogleDrive/GithubDesktopRepo/Bioinformatics"
setwd(proj_root)

if (!requireNamespace("SeuratData", quietly = TRUE)) {
  message("Installing remotes and SeuratData...")
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("satijalab/seurat-data")
}

source("scripts/seurat_generate_figures.R")
