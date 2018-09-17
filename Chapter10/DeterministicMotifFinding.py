# -*- coding: utf-8 -*-

from MySeq import MySeq


class DeterministicMotifFinding:
    """ Class for deterministic motif finding. """
    
    def __init__(self, size = 8, seqs = None):
        self.motif_size = size
        if (seqs is not None):
            self.seqs = seqs
            self.alphabet = seqs[0].alphabet()
        else:
            self.seqs = []
            self.alphabet = "ACGT" # default: DNA
        
    def __len__ (self):
        return len(self.seqs)
    
    def __getitem__(self, n):
        return self.seqs[n]
    
    def seq_size (self, i):
        return len(self.seqs[i])
    
    def read_file(self, fic, t):
        for s in open(fic, "r"):
            self.seqs.append(MySeq(s.strip().upper(),t))
        self.alphabet = self.seqs[0].alphabet()
        print(self.alphabet)
        
    def create_motif_from_indexes(self, indexes):
        #pseqs = []
        res = [[0]*self.motif_size for i in range(len(self.alphabet))]
        for i,ind in enumerate(indexes):
            subseq = self.seqs[i][ind:(ind+self.motif_size)]
            for i in range(self.motif_size):
                for k in range(len(self.alphabet)):
                    if subseq[i] == self.alphabet[k]:
                        res[k][i] = res[k][i] + 1
        return res    

    def score(self, s):
        score = 0
        mat = self.create_motif_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score += maxcol
        return score
   
    def score_multiplicative(self, s):
        score = 1.0
        mat = self.create_motif_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for  i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score *= maxcol
        return score

    # EXHAUSTIVE SEARCH 
    def next_solution (self, s):
        next_sol= [0]*len(s)
        pos = len(s) - 1     
        while pos >=0 and s[pos] == self.seq_size(pos) - self.motif_size:
            pos -= 1
        if (pos < 0): 
            next_sol = None
        else:
            for i in range(pos): 
                next_sol[i] = s[i]
            next_sol[pos] = s[pos]+1;
            for i in range(pos+1, len(s)):
                next_sol[i] = 0
        return next_sol

    def exhaustive_search(self):
        best_score = -1
        res = []
        s = [0]* len(self.seqs)
        while (s!= None):
            sc = self.score(s)
            if (sc > best_score):
                best_score = sc
                res = s
            s = self.next_solution(s)
        return res


    # # BRANCH AND BOUND     
    def next_vertex (self, s):
        res =  []
        if len(s) < len(self.seqs): # internal node -> down one level
            for i in range(len(s)): 
                res.append(s[i])
            res.append(0)
        else: # bypass
            pos = len(s)-1 
            while pos >=0 and s[pos] == self.seq_size(pos) - self.motif_size:
                pos -= 1
            if pos < 0: res = None # last solution
            else:
                for i in range(pos): res.append(s[i])
                res.append(s[pos]+1)
        return res
    
    
    def bypass (self, s):
        res =  []
        pos = len(s)-1
        while pos >=0 and s[pos] == self.seq_size(pos) - self.motif_size:
            pos -= 1
        if pos < 0: res = None 
        else:
            for i in range(pos): res.append(s[i])
            res.append(s[pos]+1)
        return res
     
    def branch_and_bound (self):
        best_score = -1
        best_motif = None
        size = len(self.seqs)
        s = [0]*size
        while s != None:
            if len(s) < size:
                # estimate the bound for current internal node
                # test if the best score can be reached
                optimum_score = self.score(s) + (size-len(s)) * self.motif_size
                if optimum_score < best_score: s = self.bypass(s)
                else: s = self.next_vertex(s)
            else:
                # test if current leaf is a better solution
                sc = self.score(s)
                if sc > best_score:
                    best_score = sc
                    best_motif = s
                s = self.next_vertex(s)
        return best_motif


    # Consensus (heuristic)   
    def heuristic_consensus(self):
        res = [0]* len(self.seqs) 
        max_score = -1;
        partial = [0,0]
        for i in range(self.seq_size(0)-self.motif_size):
            for j in range(self.seq_size(1)-self.motif_size):
                partial[0] = i
                partial[1] = j
                sc = self.score(partial);
                if(sc > max_score):
                    max_score = sc
                    res[0] = i
                    res[1] = j
        for k in range(2, len(self.seqs)):
            partial = [0]*(k+1)
            for j in range(k):
                partial[j] = res[j]
            max_score = -1
            for i in range(self.seq_size(k)-self.motif_size):
                partial[k] = i
                sc = self.score(partial)
                if(sc > max_score):
                    max_score = sc
                    res[k] = i
        return res    


    
def test1():  
    sm = DeterministicMotifFinding()
    print(sm.alphabet)
    print (len(sm.alphabet))
    sm.read_file("exampleMotifs.txt","DNA")
    sol = [25,20,2,55,59]
    print (len(sm.alphabet))
    si = sm.create_motif_from_indexes(sol)
    print (si)  
    sa = sm.score(sol)
    print(sa)
    scm = sm.score_multiplicative(sol)
    print(scm)

    
def test2():
    seq1 = MySeq("ATAGAGCTGA","DNA")
    seq2 = MySeq("ACGTAGATGA","DNA")
    seq3 = MySeq("AAGATAGGGG","DNA")
    mf = DeterministicMotifFinding(3, [seq1,seq2,seq3])
    
    print ("Exhaustive:")
    sol = mf.exhaustive_search()
    print ("Solution: " , sol)
    print ("Score: ", mf.score(sol))
    
    print ("\nBranch and Bound:")
    sol2 = mf.branch_and_bound()
    print ("Solution: " , sol2)
    print ("Score:" , mf.score(sol2))

    print ("\nHeuristic consensus: ")
    sol3 = mf.heuristic_consensus()
    print ("Solution: " , sol3)
    print ("Score:" , mf.score(sol3))

def test3():
    mf = DeterministicMotifFinding()
    mf.read_file("exampleMotifs.txt","DNA")
    print ("Branch and Bound:")
    sol = mf.branch_and_bound()
    print ("Solution: ", sol)
    print ("Score:", mf.score(sol))
    
    print ("\nHeuristic consensus: ")
    sol2 = mf.heuristic_consensus()
    print ("Solution: " , sol2)
    print ("Score:" , mf.score(sol2)) 

test1()
print()
test2()
print()
test3()