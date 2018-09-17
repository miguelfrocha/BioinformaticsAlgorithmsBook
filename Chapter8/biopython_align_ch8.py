

from Bio import Alphabet
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


seq1 = "MHQAIFIYQIGYPLKSGYIQSIRSPEYDNW"
seq2 = "MH--IFIYQIGYALKSGYIQSIRSPEY-NW"
seq3 = "MHQAIFI-QIGYALKSGY-QSIRSPEYDNW"

seqr1 = SeqRecord(Seq(seq1,Alphabet.Gapped(IUPAC.protein)),id="seq1")
seqr2 = SeqRecord(Seq(seq2,Alphabet.Gapped(IUPAC.protein)),id="seq2")
seqr3 = SeqRecord(Seq(seq3,Alphabet.Gapped(IUPAC.protein)),id="seq3")

alin = MultipleSeqAlignment([seqr1, seqr2, seqr3])
print(alin)

print(alin[1]) # 2nd sequence
print(alin[:,2]) # 3rd column
print(alin[:,3:7])  # 4th to 7th columns (all sequences)
print(alin[0].seq[:3]) # first 3 columns of seq1
print(alin[1:3,5:12]) # sequences 2 and 3; 4th to 10th column

from Bio import AlignIO
alin2 = AlignIO.read("PF05371_seed.aln", "clustal")

print("Size:", alin2.get_alignment_length())
for record in alin2:
    print (record.seq, record.id)

AlignIO.write(alin2, "example_alin.fasta", "fasta")

AlignIO.convert("PF05371_seed.aln", "clustal", "example_alin.fasta", "fasta")

