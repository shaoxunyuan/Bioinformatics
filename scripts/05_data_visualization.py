import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns

data = {
    'Sample': ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
    'Gene1': [10.5, 11.2, 15.3, 16.1, 9.8, 10.5],
    'Gene2': [8.2, 7.9, 14.2, 13.8, 8.0, 8.2],
    'Gene3': [12.3, 12.8, 18.5, 19.2, 12.0, 12.3],
    'Gene4': [9.8, 10.1, 16.1, 15.8, 9.5, 9.8],
    'Gene5': [11.2, 11.5, 17.2, 17.8, 11.0, 11.2]
}

df = pd.DataFrame(data)

print("原始数据:")
print(df)

numeric_data = df[['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5']]
sample_names = df['Sample'].values

scaler = StandardScaler()
scaled_data = scaler.fit_transform(numeric_data)

print("\n标准化后的数据:")
print(scaled_data)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_data

print("\nPCA结果:")
print(f"解释方差比: {pca.explained_variance_ratio_}")
print(f"累积解释方差: {pca.explained_variance_ratio_.cumsum()}")

pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['Sample'] = sample_names

print("\nPCA坐标:")
print(pca_df)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

axes[0, 0].scatter(pca_df['PC1'], pca_df['PC2'], s=100, alpha=0.7)
for i, txt in enumerate(pca_df['Sample']):
    axes[0, 0].annotate(txt, (pca_df['PC1'][i], pca_df['PC2'][i]), 
                     xytext=(5, 5), textcoords='offset points')
axes[0, 0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)')
axes[0, 0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)')
axes[0, 0].set_title('PCA分析')
axes[0, 0].grid(True, alpha=0.3)

linkage_matrix = linkage(scaled_data, method='ward')
dendrogram(linkage_matrix, labels=sample_names, ax=axes[0, 1], orientation='left')
axes[0, 1].set_title('层次聚类树状图')
axes[0, 1].set_xlabel('距离')

axes[1, 0].imshow(scaled_data, cmap='viridis', aspect='auto')
axes[1, 0].set_xticks(range(len(numeric_data.columns)))
axes[1, 0].set_xticklabels(numeric_data.columns, rotation=45)
axes[1, 0].set_yticks(range(len(sample_names)))
axes[1, 0].set_yticklabels(sample_names)
axes[1, 0].set_title('热图')
axes[1, 0].set_xlabel('基因')
axes[1, 0].set_ylabel('样本')

correlation_matrix = np.corrcoef(scaled_data.T)
im = axes[1, 1].imshow(correlation_matrix, cmap='coolwarm', vmin=-1, vmax=1)
axes[1, 1].set_xticks(range(len(numeric_data.columns)))
axes[1, 1].set_xticklabels(numeric_data.columns, rotation=45)
axes[1, 1].set_yticks(range(len(numeric_data.columns)))
axes[1, 1].set_yticklabels(numeric_data.columns)
axes[1, 1].set_title('基因相关性矩阵')
plt.colorbar(im, ax=axes[1, 1])

plt.tight_layout()
out_dir = os.path.join(os.path.dirname(__file__), '..', 'docs', 'images')
os.makedirs(out_dir, exist_ok=True)
plt.savefig(os.path.join(out_dir, 'data_visualization.png'), dpi=300, bbox_inches='tight')
print("\n图表已保存到 docs/images/data_visualization.png")

print("\n基因相关性矩阵:")
print(pd.DataFrame(correlation_matrix, 
                    index=numeric_data.columns, 
                    columns=numeric_data.columns).round(3))
