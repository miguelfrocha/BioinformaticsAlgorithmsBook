# -*- coding: utf-8 -*-


from MySeq import MySeq
from MyMotifs import MyMotifs

class MotifFinding:
    
    def __init__(self, size = 8, seqs = None):
        self.motif_size = size
        if (seqs != None):
            self.seqs = seqs
            self.alphabet = seqs[0].alphabet()
        else:
            self.seqs = []
                    
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
        
        
    def create_motif_from_indexes(self, indexes):
        pseqs = []
        for i,ind in enumerate(indexes):
            pseqs.append( MySeq(self.seqs[i][ind:(ind+self.motif_size)], self.seqs[i].get_seq_biotype()) )
        return MyMotifs(pseqs)
        
        
    # SCORES (include in deterministic motif finding)
        
    def score(self, s):
        score = 0
        motif = self.create_motif_from_indexes(s)
        motif.do_counts()
        mat = motif.counts
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for  i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score += maxcol
        return score
   
    def score_multiplicative(self, s):
        score = 1.0
        motif = self.create_motif_from_indexes(s)
        motif.create_pwm()
        mat = motif.pwm
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for  i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score *= maxcol
        return score     
       
    # Probabilistic Motif Finding  
    # heuristic stochastic

    def heuristic_stochastic (self):
        from random import randint
        s = [0]* len(self.seqs) 
        for k in range(len(s)):
            s[k] = randint(0, self.seq_size(k)- self.motif_size)
        motif = self.create_motif_from_indexes(s)
        motif.create_pwm()
        sc = self.score_multiplicative(s)
        bestsol = s
        improve = True
        while(improve):
            for k in range(len(s)):
                s[k] = motif.most_probable_sequence(self.seqs[k])
            if self.score_multiplicative(s) > sc: 
                sc = self.score_multiplicative(s)
                bestsol = s
                motif = self.create_motif_from_indexes(s)
                motif.create_pwm()
            else: improve = False    
        return bestsol


    # gibbs sampling
        
    def gibbs (self, iterations = 100):
        from random import randint
        s = []
        for i in range(len(self.seqs)):
            s.append(randint(0, len(self.seqs[0]) - self.motif_size - 1))
        best_s = list(s)
        best_score = self.score_multiplicative(s)
        for it in range(iterations):
            # randomly pick a sequence
            seq_idx = randint(0, len(self.seqs)-1) 
            seq_sel = self.seqs[seq_idx]
            s.pop(seq_idx)
            removed = self.seqs.pop(seq_idx)
            motif = self.create_motif_from_indexes(s)            
            motif.create_pwm()
            self.seqs.insert(seq_idx, removed)
            r = motif.probability_all_positions(seq_sel)
            pos = self.roulette(r)
            s.insert(seq_idx, pos)
            score = self.score_multiplicative(s)
            if score > best_score:
                best_score = score
                best_s = list(s)
        return best_s
    
    def roulette(self, f):
        from random import random
        tot = 0.0
        for x in f: tot += (0.01 + x) 
        val = random() * tot
        acum = 0.0
        idx = 0
        while acum < val:
            acum += (f[idx] + 0.01)
            idx += 1
        return idx-1 

      
# tests
      
def test1():  
    sm = MotifFinding()
    sm.read_file("exampleMotifs.txt","dna")
    sol = [25,20,2,55,59]
    sa = sm.score(sol)
    print(sa)
    scm = sm.score_multiplicative(sol)
    print(scm)

def test2():
    mf = MotifFinding()
    mf.read_file("exampleMotifs.txt","dna")
    print("Heuristic stochastic")
    sol = mf.heuristic_stochastic()
    print ("Solution: " , sol)
    print ("Score:" , mf.score(sol))
    print ("Score mult:" , mf.score_multiplicative(sol))
    print("Consensus:", mf.create_motif_from_indexes(sol).consensus())
    
    sol2 = mf.gibbs(10000)
    print ("Score:" , mf.score(sol2))
    print ("Score mult:" , mf.score_multiplicative(sol2))
    print("Consensus:", mf.create_motif_from_indexes(sol2).consensus())

test1()
print()
test2()