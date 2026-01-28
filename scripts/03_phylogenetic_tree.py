from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

seq1 = SeqRecord(Seq("ATGCGATCGATCGATCG"), id="Human", description="Human sequence")
seq2 = SeqRecord(Seq("ATGCGATCGATCGATCG"), id="Mouse", description="Mouse sequence")
seq3 = SeqRecord(Seq("ATGCGATCGATCGATCG"), id="Rat", description="Rat sequence")
seq4 = SeqRecord(Seq("ATGCGATCGATCGATCG"), id="Dog", description="Dog sequence")

alignment = MultipleSeqAlignment([seq1, seq2, seq3, seq4])

print("多序列比对:")
print(alignment)

calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

print("\n距离矩阵:")
print(distance_matrix)

constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

print("\n系统发育树:")
print(tree)

print("\n树的分支:")
for clade in tree.find_clades():
    if clade.name:
        print(f"  {clade.name}: 分支长度 = {clade.branch_length}")

print("\n树的Newick格式:")
print(tree.format('newick'))
