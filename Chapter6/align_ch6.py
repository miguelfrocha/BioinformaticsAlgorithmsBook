# -*- coding: utf-8 -*-

def create_submat (match, mismatch, alphabet):
    sm = {}
    for c1 in alphabet:
        for c2 in alphabet:
            if (c1 == c2):
                sm[c1+c2] = match
            else:
                sm[c1+c2] = mismatch
    return sm

# read substitution matrix from file
def read_submat_file (filename):
    sm = {}
    f = open(filename, "r")
    line = f.readline()
    tokens = line.split("\t")
    ns = len(tokens)
    alphabet = []
    for i in range(0, ns): 
        alphabet.append(tokens[i][0])
    for i in range(0,ns):
        line = f.readline();
        tokens = line.split("\t");
        for j in range(0, len(tokens)):
            k = alphabet[i]+alphabet[j]
            sm[k] = int(tokens[j])
    return sm

# score of a position (column)
def score_pos (c1, c2, sm, g):
    if c1 == "-" or c2=="-":
        return g
    else:
        return sm[c1+c2]

# score of the whole alignment
def score_align (seq1, seq2, sm, g):
    res = 0;
    for i in range(len(seq1)):
        res += score_pos (seq1[i], seq2[i], sm, g)
    return res

def score_affinegap (seq1, seq2, sm, g, r):
    res = 0
    ingap1 = False
    ingap2 = False
    for i in range(len(seq1)):
        if seq1[i]=="-":
            if ingap1: res += r
            else:
                ingap1 = True
                res += g
        elif seq2[i]=="-":
            if ingap2: res += r
            else:
                ingap2 = True
                res += g 
        else:
            if ingap1: ingap1 = False
            if ingap2: ingap2 = False
            res += sm[seq1[i]+seq2[i]]
    return res

## global alignment 

def needleman_Wunsch (seq1, seq2, sm, g):
    S = [[0]]
    T = [[0]]
    for j in range(1, len(seq2)+1):
        S[0].append(g * j)
        T[0].append(3)
    for i in range(1, len(seq1)+1):
        S.append([g * i])
        T.append([2])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos (seq1[i], seq2[j], sm, g); 
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            S[i+1].append(max(s1, s2, s3))
            T[i+1].append(max3t(s1, s2, s3))
    return (S, T)

def max3t (v1, v2, v3):
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3

def recover_align (T, seq1, seq2):
    res = ["", ""]
    i = len(seq1)
    j = len(seq2)
    while i>0 or j>0:
        if T[i][j]==1:
            res[0] = seq1[i-1] + res[0]
            res[1] = seq2[j-1] + res[1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "-" + res[0]
            res[1] = seq2[j-1] + res[1]
            j -= 1
        else:
            res[0] = seq1[i-1] + res[0]
            res[1] = "-" + res[1]
            i -= 1
    return res

## local alignment

def smith_Waterman (seq1, seq2, sm, g):
    S = [[0]]
    T = [[0]]
    maxscore = 0
    for j in range(1, len(seq2)+1):
        S[0].append(0)
        T[0].append(0)
    for i in range(1, len(seq1)+1):
        S.append([0])
        T.append([0])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos (seq1[i], seq2[j], sm, g); 
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            b = max(s1, s2, s3)
            if b <= 0:
                S[i+1].append(0)
                T[i+1].append(0)
            else:
                S[i+1].append(b)
                T[i+1].append(max3t(s1, s2, s3))
                if b > maxscore: 
                    maxscore = b
    return (S, T, maxscore)

def recover_align_local (S, T, seq1, seq2):
    res = ["", ""]
    i, j = max_mat(S)
    while T[i][j]>0:
        if T[i][j]==1:
            res[0] = seq1[i-1] + res[0]
            res[1] = seq2[j-1] + res[1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "-" + res[0];
            res[1] = seq2[j-1] + res[1] 
            j -= 1
        elif T[i][j] == 2:
            res[0] = seq1[i-1] + res[0]
            res[1] = "-" + res[1]
            i -= 1
    return res

def max_mat(mat):
    maxval = mat[0][0]
    maxrow = 0
    maxcol = 0
    for i in range(0,len(mat)):
        for j in range(0, len(mat[i])):
            if mat[i][j] > maxval:
                maxval = mat[i][j]
                maxrow = i
                maxcol = j
    return (maxrow,maxcol)

def identity(seq1, seq2, alphabet = "ACGT"):
    sm = create_submat(1,0,alphabet) 
    S,_ = needleman_Wunsch(seq1, seq2, sm, 0)
    equal = S[len(seq1)][len(seq2)]
    return equal / max(len(seq1), len(seq2))  

def edit_distance(seq1, seq2, alphabet = "ACTG"):
    sm = create_submat(0, -1, alphabet) 
    S = needleman_Wunsch(seq1, seq2,sm,-1)[0]
    res = -1*S[len(seq1)][len(seq2)]
    return res

def longest_common_subseq (seq1, seq2, alphabet = "ACGT"):
    sm = create_submat(1, 0, alphabet) 
    _,T = needleman_Wunsch(seq1, seq2, sm, 0)
    alin = recover_align(T, seq1, seq2)
    
    sizeal = len(alin[0])
    lcs = ""
    for i in range(sizeal):
        if alin[0][i] == alin[1][i]:
            lcs += alin[0][i]
    return lcs

def longest_common_string (seq1, seq2, alphabet = "ACGT"):
    m = max(len(seq1), len(seq2))
    pen = -1 * (m+1)
    sm = create_submat(1, pen, alphabet)
    S,T,_ = smith_Waterman(seq1, seq2, sm, pen)
    alinL= recover_align_local(S, T, seq1, seq2)
    return alinL[0]

def print_mat (mat):
    for i in range(0, len(mat)):
        print(mat[i]) 

## chapter 7
def align_query (query, ls, ms, g):
    bestScore = -1
    bestSeq = None
    bestAl = None
    for seq in ls:
        al = smith_Waterman(query, seq, ms, g)
        if al[2] > bestScore:
            bestScore = al[2]
            bestSeq = seq
            bestAl = al
    bestAlin = recover_align_local(bestAl[0], bestAl[1], query, bestSeq)
    return bestAlin, bestScore


### tests

def test_DNA():
    sm = create_submat(2,-2,"ACGT")
    seq1 = "-CAGTGCATG-ACATA"
    seq2 = "TCAG-GC-TCTACAGA"
    g = -3
    print(score_align(seq1, seq2, sm, g))

def test_prot():
    sm = read_submat_file("blosum62.mat")
    seq1 = "LGPSSGCASRIWTKSA"
    seq2 = "TGPS-G--S-IWSKSG"
    g = -8
    r = -2
    print(score_align(seq1, seq2, sm, g))
    print(score_affinegap(seq1, seq2, sm, g, r))
  
def test_global_alig():
    sm = read_submat_file("blosum62.mat")
    seq1 = "PHSWG"
    seq2 = "HGWAG"
    res = needleman_Wunsch(seq1, seq2, sm, -8)
    S = res[0]
    T = res[1]
    print("Score of optimal alignment:", S[len(seq1)][len(seq2)])
    print_mat(S)
    print_mat(T)
    alig = recover_align(T, seq1, seq2)
    print(alig[0])
    print(alig[1])

def test_local_alig():
    sm = read_submat_file("blosum62.mat")
    seq1 = "PHSWG"
    seq2 = "HGWAG"
    res = smith_Waterman(seq1, seq2, sm, -8)
    S = res[0]
    T = res[1]
    print("Score of optimal alignment:", res[2])
    print_mat(S)
    print_mat(T)
    alinL= recover_align_local(S, T, seq1, seq2)
    print(alinL[0])
    print(alinL[1])
 
test_DNA()
test_prot()
test_global_alig()
test_local_alig()
print(longest_common_string("ATAGAT","TTAGT"))