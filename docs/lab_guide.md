# 生物信息学实验指导

## 实验课程概述

本实验课程旨在通过实际操作，帮助学生掌握生物信息学的基本概念、方法和工具。课程涵盖序列分析、基因表达分析、系统发育分析、数据可视化和基因功能注释等核心内容。

## 实验环境准备

### 1. 软件安装

#### Python环境
```bash
# 安装Anaconda或Miniconda
# 下载地址: https://www.anaconda.com/products/distribution

# 创建虚拟环境
conda create -n bioinfo python=3.9

# 激活环境
conda activate bioinfo
```

#### 必要的Python包
```bash
# 安装生物信息学相关包
pip install biopython pandas numpy matplotlib seaborn scipy scikit-learn
pip install jupyter notebook plotly networkx
```

#### Jupyter Notebook
```bash
# 安装Jupyter Notebook
pip install jupyter

# 启动Jupyter Notebook
jupyter notebook
```

### 2. 数据准备

从课程平台下载以下数据文件：
- 基因表达数据
- 序列数据
- 变异数据

将数据文件放置在项目的`data/`目录下。

## 实验安排

### 实验1: 序列比对
**时间**: 第1-2周  
**目标**: 掌握序列比对的基本方法  
**内容**:
- 全局序列比对
- 局部序列比对
- 多序列比对  
**参考资料**: [assignment_1.md](assignments/assignment_1.md)

### 实验2: 基因表达分析
**时间**: 第3-4周  
**目标**: 学习基因表达数据的处理和分析  
**内容**:
- 数据预处理
- 差异表达分析
- 数据可视化  
**参考资料**: [assignment_2.md](assignments/assignment_2.md)

### 实验3: 系统发育树构建
**时间**: 第5-6周  
**目标**: 掌握系统发育分析方法  
**内容**:
- 多序列比对
- 距离矩阵计算
- 系统发育树构建  
**参考资料**: [assignment_3.md](assignments/assignment_3.md)

### 实验4: 数据可视化
**时间**: 第7-8周  
**目标**: 学习生物信息学数据的可视化方法  
**内容**:
- 基因表达数据可视化
- 多维数据降维
- 网络图可视化  
**参考资料**: [assignment_4.md](assignments/assignment_4.md)

### 实验5: 基因富集分析
**时间**: 第9-10周  
**目标**: 掌握基因功能富集分析方法  
**内容**:
- GO富集分析
- KEGG通路富集分析
- 富集结果可视化  
**参考资料**: [assignment_5.md](assignments/assignment_5.md)

### 实验6: 单细胞数据分析（Seurat）
**时间**: 第11-12周  
**目标**: 掌握单细胞 RNA 测序数据的质控、降维、聚类与可视化  
**内容**:
- Seurat 对象与 PBMC 示例数据
- 质控、归一化、高变基因与降维（PCA、UMAP）
- 聚类与细胞类型注释
- 常用图：DimPlot、FeaturePlot、VlnPlot、DotPlot、DoHeatmap 等  
**参考资料**: [07-seurat.html](07-seurat.html)

## 实验报告要求

### 报告结构
1. **实验目的**: 简述实验的目标和意义
2. **实验原理**: 介绍实验涉及的理论知识
3. **实验内容**: 详细描述实验步骤
4. **实验结果**: 展示实验结果和图表
5. **结果分析**: 对结果进行解释和讨论
6. **问题与讨论**: 遇到的问题及解决方案
7. **结论**: 总结实验收获

### 报告格式
- 使用提供的实验报告模板
- 字体：宋体，小四号
- 行距：1.5倍
- 图表：清晰标注，附图注
- 代码：关键代码片段附在报告中

### 评分标准
- 实验完成度 (40%)
- 结果准确性 (30%)
- 报告质量 (20%)
- 创新性 (10%)

## 常见问题解答

### Q1: 如何处理安装包时出现的错误？
A: 
1. 检查网络连接
2. 尝试使用国内镜像源：`pip install -i https://pypi.tuna.tsinghua.edu.cn/simple 包名`
3. 更新pip：`pip install --upgrade pip`

### Q2: Jupyter Notebook无法启动怎么办？
A:
1. 检查是否正确激活虚拟环境
2. 尝试重新安装Jupyter：`pip install --upgrade jupyter`
3. 检查端口是否被占用

### Q3: 数据文件读取失败？
A:
1. 检查文件路径是否正确
2. 确认文件格式是否正确
3. 检查文件编码（通常为UTF-8）

### Q4: 内存不足如何处理？
A:
1. 使用数据分块处理
2. 减少数据规模进行测试
3. 使用更高效的数据类型

## 学习资源

### 在线资源
- [Biopython官方文档](https://biopython.org/)
- [Python生物信息学教程](https://www.biostars.org/)
- [NCBI数据库](https://www.ncbi.nlm.nih.gov/)
- [KEGG数据库](https://www.genome.jp/kegg/)

### 推荐书籍
- 《生物信息学算法导论》
- 《Python生物信息学数据分析》
- 《生物信息学：方法与实践》

## 联系方式

如有问题，请联系：
- 助教邮箱：bioinfo@example.com
- 课程论坛：[课程论坛链接]
- 办公时间：每周三下午2-4点

## 注意事项

1. **学术诚信**: 严禁抄袭他人实验报告，引用他人成果请注明出处
2. **数据安全**: 定期备份实验数据和代码
3. **时间管理**: 合理安排实验时间，避免最后时刻匆忙完成
4. **代码规范**: 保持代码整洁，添加必要的注释
5. **实验记录**: 详细记录实验过程和遇到的问题

祝实验顺利！