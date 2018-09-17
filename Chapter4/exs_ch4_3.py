# -*- coding: utf-8 -*-

def gc_duplets(base_seq):
	seq = base_seq.upper()
	return seq.count('GC')

seq1 = 'AAGCGCTTGCG'
seq2 = 'ATCGCCGGCG'
print (gc_duplets(seq1))
print (gc_duplets(seq2))
