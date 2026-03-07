import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

genes = {
    'Gene_ID': ['GENE001', 'GENE002', 'GENE003', 'GENE004', 'GENE005', 
                'GENE006', 'GENE007', 'GENE008', 'GENE009', 'GENE010'],
    'Gene_Name': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS', 
                   'PTEN', 'AKT1', 'BCL2', 'CDK4', 'RB1'],
    'Pathway': ['DNA Repair', 'Cell Cycle', 'Growth Factor', 'Transcription', 'Signal Transduction',
                'Signal Transduction', 'Signal Transduction', 'Apoptosis', 'Cell Cycle', 'Cell Cycle'],
    'Log2FC': [2.5, 3.2, -1.8, 1.5, 2.8, -2.1, 1.9, -1.5, 2.2, -1.2],
    'P_Value': [0.001, 0.0001, 0.05, 0.01, 0.0005, 0.02, 0.008, 0.03, 0.002, 0.04]
}

df = pd.DataFrame(genes)

print("基因表达数据:")
print(df)

df['-log10_P'] = -np.log10(df['P_Value'])

print("\n基因富集分析:")

pathway_counts = df.groupby('Pathway').size()
print("\n各通路基因数:")
print(pathway_counts)

pathway_stats = df.groupby('Pathway').agg({
    'Log2FC': ['mean', 'std'],
    'P_Value': 'min'
}).round(4)

print("\n各通路统计:")
print(pathway_stats")

upregulated = df[df['Log2FC"] > 1]
downregulated = df[df['Log2FC"] < -1]

print(f"\n上调基因: {len(upregulated)}个")
print(upregulated[['Gene_ID", "Gene_Name", "Log2FC", "P_Value"]])

print(f"\n下调基因: {len(downregulated)}个")
print(downregulated[['Gene_ID", "Gene_Name", "Log2FC", "P_Value"]])

fig, axes = plt.subplots(2, 2, figsize=(14, 12"))

axes[0, 0].scatter(df['Log2FC"], df['-log10_P"], s=100, alpha=0.7)
for i, txt in enumerate(df['Gene_Name']):
    axes[0, 0].annotate(txt, (df['Log2FC'][i], df['-log10_P'][i]), 
                       xytext=(5, 5), textcoords='offset points', fontsize=8)
axes[0, 0].axhline(y=-np.log10(0.05), color='r', linestyle='--', label='p=0.05')
axes[0, 0].axvline(x=1, color='r', linestyle='--', label='FC=2')
axes[0, 0].axvline(x=-1, color='r', linestyle='--')
axes[0, 0].set_xlabel('Log2 Fold Change')
axes[0, 0].set_ylabel('-log10(p-value)')
axes[0, 0].set_title('火山图')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

pathway_counts.plot(kind='bar', ax=axes[0, 1], color='steelblue', alpha=0.8)
axes[0, 1].set_xlabel('通路')
axes[0, 1].set_ylabel('基因数')
axes[0, 1].set_title('各通路基因数')
axes[0, 1].tick_params(axis='x', rotation=45)
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].bar(df['Gene_Name'], df['Log2FC'], color=['green' if x > 0 else 'red' for x in df['Log2FC']], alpha=0.7)
axes[1, 0].axhline(y=1, color='r', linestyle='--', label='FC=2')
axes[1, 0].axhline(y=-1, color='r', linestyle='--')
axes[1, 0].set_xlabel('基因')
axes[1, 0].set_ylabel('Log2 Fold Change')
axes[1, 0].set_title('基因表达变化')
axes[1, 0].tick_params(axis='x', rotation=45)
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].bar(df['Gene_Name'], -np.log10(df['P_Value']), color='orange', alpha=0.7)
axes[1, 1].axhline(y=-np.log10(0.05), color='r', linestyle='--', label='p=0.05')
axes[1, 1].set_xlabel('基因')
axes[1, 1].set_ylabel('-log10(p-value)')
axes[1, 1].set_title('基因显著性')
axes[1, 1].tick_params(axis='x', rotation=45)
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
out_dir = os.path.join(os.path.dirname(__file__), '..', 'docs', 'images')
os.makedirs(out_dir, exist_ok=True)
plt.savefig(os.path.join(out_dir, 'gene_enrichment.png'), dpi=300, bbox_inches='tight')
print("\n图表已保存到 docs/images/gene_enrichment.png")

print("\n富集分析结果:")
print(f"显著差异表达基因数: {len(df[df['P_Value'] < 0.05])}")
print(f"上调基因数: {len(upregulated)}")
print(f"下调基因数: {len(downregulated)}")

print("\n各通路详细统计:")
for pathway in df['Pathway'].unique():
    pathway_genes = df[df['Pathway'] == pathway]
    print(f"\n{pathway}:")
    print(f"  基因数: {len(pathway_genes)}")
    print(f"  平均Log2FC: {pathway_genes['Log2FC'].mean():.2f}")
    print(f"  最小P值: {pathway_genes['P_Value'].min():.4f}")
    print(f"  基因列表: {', '.join(pathway_genes['Gene_Name'].tolist())}")
