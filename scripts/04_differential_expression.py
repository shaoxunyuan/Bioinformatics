import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

data = {
    'Gene': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5', 'Gene6', 'Gene7', 'Gene8', 'Gene9', 'Gene10'],
    'Control_1': [10.5, 8.2, 12.3, 9.8, 11.2, 10.1, 9.5, 10.8, 11.0, 9.9],
    'Control_2': [11.2, 7.9, 12.8, 10.1, 11.5, 10.3, 9.8, 11.0, 11.3, 10.2],
    'Treatment_1': [15.3, 14.2, 18.5, 16.1, 17.2, 15.8, 14.9, 16.5, 17.0, 15.5],
    'Treatment_2': [16.1, 13.8, 19.2, 15.8, 17.8, 16.2, 15.2, 16.8, 17.5, 15.9]
}

df = pd.DataFrame(data)

print("基因表达数据:")
print(df)

control_mean = df[['Control_1', 'Control_2']].mean(axis=1)
treatment_mean = df[['Treatment_1', 'Treatment_2']].mean(axis=1)

df['Control_Mean'] = control_mean
df['Treatment_Mean'] = treatment_mean
df['Fold_Change'] = treatment_mean / control_mean
df['Log2_Fold_Change'] = np.log2(df['Fold_Change'])

print("\n统计分析:")
print(df[['Gene', 'Control_Mean', 'Treatment_Mean', 'Fold_Change', 'Log2_Fold_Change']])

control_std = df[['Control_1', 'Control_2']].std(axis=1)
treatment_std = df[['Treatment_1', 'Treatment_2']].std(axis=1)

n1 = 2
n2 = 2

pooled_std = np.sqrt(((n1-1)*control_std**2 + (n2-1)*treatment_std**2) / (n1+n2-2))
se = pooled_std * np.sqrt(1/n1 + 1/n2)

t_stat = (treatment_mean - control_mean) / se
df_degrees = n1 + n2 - 2
p_values = 2 * (1 - stats.t.cdf(np.abs(t_stat), df_degrees))

df['t_statistic'] = t_stat
df['p_value'] = p_values
df['-log10(p_value)'] = -np.log10(p_values)

print("\n差异表达分析:")
print(df[['Gene', 't_statistic', 'p_value', '-log10(p_value)']])

upregulated = df[df['Fold_Change'] > 1.5]
downregulated = df[df['Fold_Change'] < 0.67]

print(f"\n上调基因数: {len(upregulated)}")
print(f"下调基因数: {len(downregulated)}")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

axes[0, 0].scatter(df['Log2_Fold_Change'], df['-log10(p_value)'], alpha=0.6)
axes[0, 0].axhline(y=-np.log10(0.05), color='r', linestyle='--', label='p=0.05')
axes[0, 0].axvline(x=np.log2(1.5), color='r', linestyle='--', label='FC=1.5')
axes[0, 0].axvline(x=np.log2(0.67), color='r', linestyle='--', label='FC=0.67')
axes[0, 0].set_xlabel('Log2 Fold Change')
axes[0, 0].set_ylabel('-log10(p-value)')
axes[0, 0].set_title('火山图')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

x = np.arange(len(df))
width = 0.35
axes[0, 1].bar(x - width/2, df['Control_Mean'], width, label='Control', alpha=0.8)
axes[0, 1].bar(x + width/2, df['Treatment_Mean'], width, label='Treatment', alpha=0.8)
axes[0, 1].set_xlabel('Gene')
axes[0, 1].set_ylabel('Expression Level')
axes[0, 1].set_title('基因表达水平比较')
axes[0, 1].set_xticks(x)
axes[0, 1].set_xticklabels(df['Gene'], rotation=45)
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].hist(df['Log2_Fold_Change'], bins=10, edgecolor='black', alpha=0.7)
axes[1, 0].set_xlabel('Log2 Fold Change')
axes[1, 0].set_ylabel('Frequency')
axes[1, 0].set_title('Fold Change 分布')
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].scatter(df['Control_Mean'], df['Treatment_Mean'], alpha=0.6)
max_val = max(df['Control_Mean'].max(), df['Treatment_Mean'].max())
axes[1, 1].plot([0, max_val], [0, max_val], 'r--', label='1:1 line')
axes[1, 1].set_xlabel('Control Mean')
axes[1, 1].set_ylabel('Treatment Mean')
axes[1, 1].set_title('MA图')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
out_dir = os.path.join(os.path.dirname(__file__), '..', 'docs', 'images')
os.makedirs(out_dir, exist_ok=True)
plt.savefig(os.path.join(out_dir, 'expression_analysis.png'), dpi=300, bbox_inches='tight')
print("\n图表已保存到 docs/images/expression_analysis.png")
