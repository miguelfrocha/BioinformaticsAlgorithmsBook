#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def aprox_compare(seq1, seq2, d):
    mismatches = 0
    i = 0
    while mismatches <= d and i < len(seq1): # assuming equal len of the sequences
        if seq1[i] != seq2[i]:
            mismatches += 1
        i = i + 1
    if mismatches <= d: return True
    else: return False

def aprox_match(seq, p, d):
    res = []
    for pos in range(len(seq)-len(p)+1):
        if aprox_compare( seq[pos:pos+len(p)] , p, d):
            res.append(pos)
    return res

print( aprox_compare("ATGAT","ATTAA",1))
print( aprox_compare("ATGAT","ATTAA",2))
