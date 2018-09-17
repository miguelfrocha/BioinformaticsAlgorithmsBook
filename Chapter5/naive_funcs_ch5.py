# -*- coding: utf-8 -*-

def search_first_occ(seq, pattern):
    found = False
    i = 0
    while i <= len(seq)-len(pattern) and not found:
        j = 0
        while j < len(pattern) and pattern[j]==seq[i+j]: 
            j = j + 1
        if j== len(pattern): found = True
        else: i += 1
    if found: return i
    else: return -1
    
def search_all_occurrences(seq, pattern):
    res = [] 
    for i in range(len(seq)-len(pattern)+1):
        j = 0
        while j < len(pattern) and pattern[j]==seq[i+j]:
            j = j + 1
        if j == len(pattern): 
            res.append(i)
    return res

def test_pat_search():
    seq = input("Input sequence: ")
    pat = input("Input pattern: ")
    print(pat, "occurs in the following positions:", )
    print( search_all_occurrences(seq, pat) )
    
test_pat_search()    

#seqDNA = "ATAGAATAGATAATAGTC"
#print( search_first_occ(seqDNA, "GAAT") )
#print( search_first_occ(seqDNA, "TATA") )
#print( search_all_occurrences(seqDNA, "AAT") )


