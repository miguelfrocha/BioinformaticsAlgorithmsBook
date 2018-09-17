# -*- coding: utf-8 -*-

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo 

alignments = pairwise2.align.globalxx("ATAGAGAATAG", "ATGGCAGATAGA")

print (len(alignments))

for a in alignments: 
    print(format_alignment(*a))

matrix = MatrixInfo.blosum62
for a in pairwise2.align.globalds("KEVLA", "EVSAW", matrix,-4,-1):
    print(format_alignment(*a))
    

local_dna = pairwise2.align.localms("ATAGAGAATAG", "GGGAGAATC", 3,-2,-3,-3)
for a in local_dna: print(format_alignment(*a))

local_prot = pairwise2.align.localds("KEVLA", "EVSAW", matrix,-4,-1)
for a in local_prot: print(format_alignment(*a))
