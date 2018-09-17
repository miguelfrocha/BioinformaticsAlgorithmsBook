#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def subseqs(seqDNA, seqAA):

    d = {"GCT":"A", "GCC":"A", "GCA":"A", "GCC":"A", 
          "TGT":"C", "TGC":"C",
          "GAT":"D", "GAC":"D",
          "GAA":"E", "GAG":"E",
          "TTT":"F", "TTC":"F",
          "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
          "CAT":"H", "CAC":"H",
          "ATA":"I", "ATT":"I", "ATC":"I",
          "AAA":"K", "AAG":"K",
          "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
          "ATG":"M", "AAT":"N", "AAC":"N",
          "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
          "CAA":"Q", "CAG":"Q",
          "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
          "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
          "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
          "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
          "TGG":"W",
          "TAT":"Y", "TAC":"Y",
          "TAA":"_", "TAG":"_", "TGA":"_"
          }
    subseqs = []
    lastpos = len(seqDNA) - len(seqAA)*3 + 1
    for i in range(lastpos): # initial position
        prot = ""
        for k in range(i, i + len(seqAA)*3, 3):
            codao = seqDNA[k : k+3]
            aa = d[codao]
            prot += aa
        if prot == seqAA:
            subseqs.append(i)
    return subseqs

## main 
seqAA = input("Insert aminoacid sequence:")
seqDNA = input("Insert DNA sequence:")
subs = subseqs(seqDNA, seqAA)
if subs:
    print("Subsequences encoding the target aminoacid sequence")
    for s in subs:
        print(seqDNA[s:s+len(seqAA*3)])
        
else:
    print("No subsequences encode the target aminoacid sequence")


