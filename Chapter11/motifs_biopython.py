# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC

instances = []
instances.append(Seq("TATAA",IUPAC.unambiguous_dna))
instances.append(Seq("TATTA",IUPAC.unambiguous_dna))
instances.append(Seq("TTTAT",IUPAC.unambiguous_dna))
instances.append(Seq("TATAC",IUPAC.unambiguous_dna))

m = motifs.create(instances, IUPAC.unambiguous_dna)

print(type(m))
print(m)
print(len(m))
print(m.consensus)
print(m.pwm)
print(m.counts)
print(m.pssm)

m.weblogo("mymotif.png")

pwm = m.counts.normalize(pseudocounts=0.5)
pssm = pwm.log_odds()
print(pwm)
print(pssm)


# exact matches of the instances

test_seq=Seq("TTTTATACACTGCATATAACAACCCAAGCATTATAA", IUPAC.unambiguous_dna)

for pos, seq in m.instances.search(test_seq):
    print (pos, " ", seq)
    
# using PSSM to score matches
for position, score in pssm.search(test_seq, threshold=4.0):
    print ("Position %d: score = %5.3f" % (position, score) )
    
# scores for all positions
print(pssm.calculate(test_seq))

