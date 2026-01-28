from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

dna_seq = Seq("ATGCGATCGATCGATCG")

print("DNA序列:", dna_seq)
print("序列长度:", len(dna_seq))
print("反向互补序列:", dna_seq.reverse_complement())

rna_seq = dna_seq.transcribe()
print("\n转录为RNA:", rna_seq)
print("RNA反向互补:", rna_seq.reverse_complement())

protein_seq = dna_seq.translate()
print("\n翻译为蛋白质:", protein_seq

print("\n序列组成分析:")
print("A:", dna_seq.count('A'))
print("T:", dna_seq.count('T'))
print("G:", dna_seq.count('G'))
print("C:", dna_seq.count('C'))

gc_content = (dna_seq.count('G') + dna_seq.count('C')) / len(dna_seq) * 100
print(f"GC含量: {gc_content:.2f}%")

print("\n序列操作:")
print("前5个碱基:", dna_seq[:5])
print("后5个碱基:", dna_seq[-5:])
print("查找'ATG':", dna_seq.find('ATG')))

record = SeqRecord(
    dna_seq,
    id='gene001',
    description='Example gene sequence'
)

print("\n序列记录:")
print(f"ID: {record.id}")
print(f"描述: {record.description}")
print(f"序列: {record.seq}")

print("\n序列格式转换:")
print("FASTA格式:")
print(record.format("fasta"))

print("\n序列比对准备:")
seq1 = Seq("ATGCGATCGATCGATCG")
seq2 = Seq("ATGCGATCGATCGATCG")

print(f"序列1: {seq1}")
print(f"序列2: {seq2}")
print(f"序列相等: {seq1 == seq2}")

print("\n序列变异分析:")
mutations = []
for i in range(len(seq1)):
    if seq1[i] != seq2[i]:
        mutations.append((i+1, seq1[i], seq2[i]))

if mutations:
    print("检测到变异:")
    for pos, ref, alt in mutations:
        print(f"  位置 {pos}: {ref} -> {alt}")
else:
    print("未检测到变异")

print("\n密码子分析:")
codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]
print("密码子列表:", codons)
print("密码子数:", len(codons))

start_codon = "ATG"
stop_codons = ["TAA", "TAG", "TGA"]

if codons[0] == start_codon:
    print("起始密码子: ATG (正确)")
else:
    print("起始密码子: 未找到")

if codons[-1] in stop_codons:
    print("终止密码子:", codons[-1], "(正确)")
else:
    print("终止密码子:", codons[-1] if codons else "无", "(非标准)")

print("\n序列特征:")
print("AT含量:", (dna_seq.count('A') + dna_seq.count('T')) / len(dna_seq) * 100, "%")
print("GC含量:", gc_content, "%")
print("序列复杂度:", len(set(str(dna_seq))) / len(dna_seq) * 100, "%")

print("\n序列统计:")
print("总碱基数:", len(dna_seq))
print("嘌呤数 (A+G):", dna_seq.count('A') + dna_seq.count('G'))
print("嘧啶数 (T+C):", dna_seq.count('T') + dna_seq.count('C'))
