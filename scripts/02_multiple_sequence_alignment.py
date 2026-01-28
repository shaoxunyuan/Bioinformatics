from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

seq1 = SeqRecord(Seq("ATGCGATCGATCGATCG"), id="seq1", description="Sequence 1")
seq2 = SeqRecord(Seq("ATGCGATCGATCGATCG"), id="seq2", description="Sequence 2")
seq3 = SeqRecord(Seq("ATGCGATCGATCGATCG"), id="seq3", description="Sequence 3")

alignment = MultipleSeqAlignment([seq1, seq2, seq3])

print("多序列比对:")
print(alignment)

print("\n比对长度:", alignment.get_alignment_length())

print("\n比对统计:")
print(f"序列数: {len(alignment)}")
print(f"比对长度: {alignment.get_alignment_length()}")

print("\n各序列:")
for record in alignment:
    print(f"ID: {record.id}, 长度: {len(record.seq)}, 序列: {record.seq}")

consensus = alignment.consensus
print(f"\n一致性序列: {consensus}")

gc_content = (consensus.count('G') + consensus.count('C')) / len(consensus) * 100
print(f"GC含量: {gc_content:.2f}%")
