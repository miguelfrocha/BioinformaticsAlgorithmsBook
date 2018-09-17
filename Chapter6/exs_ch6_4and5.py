#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from align_ch6 import score_pos, create_submat, needleman_Wunsch, smith_Waterman, recover_align, print_mat

def needlemanWunsch_with_ties (seq1, seq2, sm, g):
    S = [[0]]
    T = [[0]]
    for j in range(1, len(seq2)+1):
        S[0].append(g * j)
        T[0].append([3])
    for i in range(1, len(seq1)+1):
        S.append([g * i])
        T.append([[2]])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos (seq1[i], seq2[j], sm, g); 
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            S[i+1].append(max(s1, s2, s3))
            T[i+1].append(max3t_with_ties(s1, s2, s3))
    return (S, T)


def smithWaterman_with_ties(seq1, seq2, sm, g):
    S = [[0]]
    T = [[0]]
    maxscore = 0
    for j in range(1, len(seq2)+1):
        S[0].append(0)
        T[0].append([0])
    for i in range(1, len(seq1)+1):
        S.append([0])
        T.append([[0]])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos (seq1[i], seq2[j], sm, g); 
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            b = max(s1, s2, s3)
            if b <= 0:
                S[i+1].append(0)
                T[i+1].append([0])
            else:
                S[i+1].append(b)
                T[i+1].append(max3t_with_ties(s1, s2, s3))
                if b > maxscore: 
                    maxscore = b
    return (S, T, maxscore)

def max3t_with_ties(v1, v2, v3):
    if v1 > v2:
        if v1 > v3: 
            return [1]
        elif v1 == v3:
            return [1,3]
        else:
            return [3]
    elif v1 == v2:
        if v1 > v3: 
            return [1,2]
        elif v1 == v3:
            return [1,2,3]
        else:
            return [3]
    else:
        if v2 > v3: return [2]
        elif v2 == v3:
            return [2,3]
        else: return [3]

def recover_align_with_ties (T, seq1, seq2):
    i = len(seq1)
    j = len(seq2)  
    alins = [["", "", i,j]]
    res = []
    while alins:
        al = alins.pop(0)
        i = al[2]
        j = al[3]
        if i==0 and j==0:
            res.append(al[:2])
        else:
            for t in T[i][j]:
                p = []
                if t==1:
                    p.append(seq1[i-1] + al[0])
                    p.append(seq2[j-1] + al[1])
                    p.append(i-1)
                    p.append(j-1)
                elif t == 3:
                    p.append("-" + al[0])
                    p.append(seq2[j-1] + al[1])
                    p.append(i)
                    p.append(j-1)
                else:
                    p.append(seq1[i-1] + al[0])
                    p.append("-" + al[1])
                    p.append(i-1)
                    p.append(j)
                alins.append(p)
    return res

def recoverAlignLocal_with_ties (S, T, seq1, seq2):
    maxval = S[0][0]
    maxtups = []
    for i in range(0,len(S)):
        for j in range(0, len(S[i])):
            if S[i][j] > maxval:
                maxval = S[i][j]
                maxtups = [(i,j)]
            elif S[i][j] == maxval:
                maxtups.append((i,j))
    alins = []
    for (i,j) in maxtups:
        alins.append(["", "", i,j])        
    res = []
    while alins:
        al = alins.pop(0)
        i = al[2]
        j = al[3]   
        if (i==0 and j==0) or (0 in T[i][j]):
            res.append(al[:2])
        else:
            for t in T[i][j]:
                p = []
                if t==1:
                    p.append(seq1[i-1] + al[0])
                    p.append(seq2[j-1] + al[1])
                    p.append(i-1)
                    p.append(j-1)
                elif t == 3:
                    p.append("-" + al[0])
                    p.append(seq2[j-1] + al[1])
                    p.append(i)
                    p.append(j-1)
                else:
                    p.append(seq1[i-1] + al[0])
                    p.append("-" + al[1])
                    p.append(i-1)
                    p.append(j)
                alins.append(p)
    return res


def test_global():
    seq1 = "TAGAAT"
    seq2 = "TAAGAT"
    sm = create_submat(1,0,"ACGT")

    res = needleman_Wunsch(seq1, seq2, sm, 0)
    S = res[0]
    T = res[1]
    print_mat(S)
    print()
    print_mat(T)
    print()
    
    _,T = needlemanWunsch_with_ties(seq1, seq2, sm, 0)
    alins = recover_align_with_ties(T, seq1, seq2)
    
    print("Score of optimal alignment:", S[len(seq1)][len(seq2)])
    
    if len(alins)==1: print("No alternative alignments")
    print("Alternative optimal global alignments:")
    print(alins)

def test_local():
    seq1 = "ATAGAT"
    seq2 = "TAAGTC"
    sm = create_submat(1,0,"ACGT")

    res = smith_Waterman(seq1, seq2, sm, 0)
    S = res[0]
    T = res[1]
    print_mat(S)
    print()
    print_mat(T)
    print()
    
    S,T,_ = smithWaterman_with_ties(seq1, seq2, sm, 0)
    alins = recoverAlignLocal_with_ties(S, T, seq1, seq2)
    
    print("Score of optimal local alignment:", res[2])
    
    if len(alins)==1: print("No alternative alignments")
    print("Alternative optimal global alignments:")
    print(alins)
   
    
test_global()
test_local()