南京中医药大学 人工智能与信息技术学院

https://shaoxunyuan.github.io/Bioinformatics/

# 生物信息学实验课程

本仓库包含生物信息学实验课程的教学材料，包括参考文档、实验作业、Python脚本、Jupyter Notebook教程和示例数据。

## 课程内容

### 📚 参考文档
- [序列比对原理](docs/01-alignment.html)
- [系统发育分析原理](docs/02-evolution.html)
- [序列映射原理](docs/03-mapping.html)
- [基因表达分析原理](docs/04-expression.html)
- [数据可视化原理](docs/05-visualization.html)
- [基因富集分析原理](docs/06-enrichments.html)
- [实验指导文档](docs/lab_guide.md)

### 📝 实验作业
- [实验作业 1: 序列比对](assignments/assignment_1.md) - 全局和局部序列比对
- [实验作业 2: 基因表达分析](assignments/assignment_2.md) - 差异表达分析和数据可视化
- [实验作业 3: 系统发育树构建](assignments/assignment_3.md) - 多序列比对和系统发育分析
- [实验作业 4: 数据可视化](assignments/assignment_4.md) - 生物信息学数据可视化方法
- [实验作业 5: 基因富集分析](assignments/assignment_5.md) - GO和KEGG通路富集分析

详见 [作业说明](assignments/README.md)

### 💻 Python脚本
- [01_sequence_alignment.py](scripts/01_sequence_alignment.py) - 序列比对脚本
- [02_multiple_sequence_alignment.py](scripts/02_multiple_sequence_alignment.py) - 多序列比对脚本
- [03_phylogenetic_tree.py](scripts/03_phylogenetic_tree.py) - 系统发育树构建脚本
- [04_differential_expression.py](scripts/04_differential_expression.py) - 差异表达分析脚本
- [05_data_visualization.py](scripts/05_data_visualization.py) - 数据可视化脚本
- [06_gene_enrichment.py](scripts/06_gene_enrichment.py) - 基因富集分析脚本
- [07_sequence_operations.py](scripts/07_sequence_operations.py) - 序列操作脚本

### 📓 Jupyter Notebook教程
- [01_sequence_alignment.ipynb](notebooks/01_sequence_alignment.ipynb) - 序列比对交互教程
- [02_gene_expression.ipynb](notebooks/02_gene_expression.ipynb) - 基因表达分析交互教程
- [03_data_visualization.ipynb](notebooks/03_data_visualization.ipynb) - 数据可视化交互教程

### 📊 示例数据
- [sample_sequences.fasta](data/sample_sequences.fasta) - 示例DNA序列
- [sample_proteins.fasta](data/sample_proteins.fasta) - 示例蛋白质序列
- [sample_rnaseq_counts.txt](data/sample_rnaseq_counts.txt) - 示例RNA-seq计数数据
- [sample_variants.txt](data/sample_variants.txt) - 示例变异数据

## 快速开始

### 环境配置
```bash
# 创建虚拟环境
conda create -n bioinfo python=3.9
conda activate bioinfo

# 安装必要的包
pip install biopython pandas numpy matplotlib seaborn scipy scikit-learn
pip install jupyter notebook plotly networkx
```

### 运行Jupyter Notebook
```bash
jupyter notebook
```

### 运行Python脚本
```bash
python scripts/01_sequence_alignment.py
```

## 实验安排

| 周次 | 实验内容 | 作业 |
|------|---------|------|
| 1-2  | 序列比对 | [作业1](assignments/assignment_1.md) |
| 3-4  | 基因表达分析 | [作业2](assignments/assignment_2.md) |
| 5-6  | 系统发育树构建 | [作业3](assignments/assignment_3.md) |
| 7-8  | 数据可视化 | [作业4](assignments/assignment_4.md) |
| 9-10 | 基因富集分析 | [作业5](assignments/assignment_5.md) |

## 学习资源

### 在线资源
- [Biopython官方文档](https://biopython.org/)
- [NCBI数据库](https://www.ncbi.nlm.nih.gov/)
- [KEGG数据库](https://www.genome.jp/kegg/)

### 推荐书籍
- 《生物信息学算法导论》
- 《Python生物信息学数据分析》
- 《生物信息学：方法与实践》

## 联系方式

- 助教邮箱：bioinfo@example.com
- 课程论坛：[课程论坛链接]
- 办公时间：每周三下午2-4点

## 许可证

本项目采用 MIT 许可证。详见 [LICENSE](LICENSE) 文件。
