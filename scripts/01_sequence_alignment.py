from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq

seq1 = "ATGCGATCGATCGATCG"
seq2 = "ATGCGATCGATCGATCG"

alignments = pairwise2.align.globalxx(seq1, seq2)

print("全局比对结果:")
for alignment in alignments[:1]:
    print(format_alignment(*alignment))

seq3 = "ATGCGATCG"
seq4 = "ATGCGATCGATCGATCGATCGATCG"

local_alignments = pairwise2.align.localxx(seq3, seq4)

print("\n局部比对结果:")
for alignment in local_alignments[:1]:
    print(format_alignment(*alignment))

seq5 = "ATGCGATCGATCGATCG"
seq6 = "ATGCGATCGATCGATCG"

score_alignments = pairwise2.align.globalms(seq5, seq6, 2, -1, -0.5, -0.1)

print("\n带得分的全局比对结果:")
for alignment in score_alignments[:1]:
    print(format_alignment(*alignment))
