# -*- coding: utf-8 -*-


def create_matrix_zeros (nrows, ncols):
    res = [ ] 
    for i in range(0, nrows):
        res.append([0]*ncols)
    return res

def print_matrix(mat):
    for i in range(0, len(mat)): print(mat[i])

class MyMotifs:

    def __init__(self, seqs = [], pwm = [], alphabet = None):
        if seqs:
            self.size = len(seqs[0])
            self.seqs = seqs # objet from class MySeq
            self.alphabet = seqs[0].alphabet()
            self.do_counts()
            self.create_pwm()
        else:
            self.pwm = pwm
            self.size = len(pwm[0])
            self.alphabet = alphabet   
        
    def __len__ (self):
        return self.size        
        
    def do_counts(self):
        self.counts = create_matrix_zeros(len(self.alphabet), self.size)
        for s in self.seqs:
            for i in range(self.size):
                lin = self.alphabet.index(s[i])
                self.counts[lin][i] += 1
                
    def create_pwm(self):
        if self.counts == None: self.do_counts()
        self.pwm = create_matrix_zeros(len(self.alphabet), self.size)
        for i in range(len(self.alphabet)):
            for j in range(self.size):
                self.pwm[i][j] = float(self.counts[i][j]) / len(self.seqs)
                
    def consensus(self):
        """ returns the sequence motif obtained with the most frequent symbol at each position of the motif"""
        res = ""
        for j in range(self.size):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alphabet) ):
                if self.counts[i][j] > maxcol: 
                    maxcol = self.counts[i][j]
                    maxcoli = i
            res += self.alphabet[maxcoli]        
        return res

    def masked_consensus(self):
        """ returns the sequence motif obtained with the symbol that occurrs in at least 50% of the input sequences"""
        res = ""
        for j in range(self.size):
            maxcol = self.counts[0][j]
            maxcoli = 0
            for i in range(1, len(self.alphabet) ):
                if self.counts[i][j] > maxcol: 
                    maxcol = self.counts[i][j]
                    maxcoli = i
            if maxcol > len(self.seqs) / 2:
                res += self.alphabet[maxcoli]        
            else:
                res += "-"
        return res

    def probability_sequence(self, seq):
        res = 1.0
        for i in range(self.size):
            lin = self.alphabet.index(seq[i])
            res *= self.pwm[lin][i]
        return res
    
    def probability_all_positions(self, seq):
        res = []
        for k in range(len(seq)-self.size+1):
            res.append(self.probability_sequence(seq))
        return res

    def most_probable_sequence(self, seq):
        """ Returns the index of the most probable sub-sequence of the input sequence"""
        maximum = -1.0
        maxind = -1
        for k in range(len(seq)-self.size):
            p = self.probability_sequence(seq[k:k+ self.size])
            if(p > maximum):
                maximum = p
                maxind = k
        return maxind

    def create_motif(self, seqs):
        from MySeq import MySeq
        l = []
        for s in seqs: 
            ind = self.most_probable_sequence(s.seq)
            subseq = MySeq ( s[ind:(ind+self.size)], s.get_seq_biotype() )
            l.append(subseq)
        
        return MyMotifs(l) 


def test():
    from MySeq import MySeq
    seq1 = MySeq("AAAGTT")
    seq2 = MySeq("CACGTG")
    seq3 = MySeq("TTGGGT")
    seq4 = MySeq("GACCGT")
    seq5 = MySeq("AACCAT")
    seq6 = MySeq("AACCCT")
    seq7 = MySeq("AAACCT")
    seq8 = MySeq("GAACCT")
    lseqs = [seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8]
    motifs = MyMotifs(lseqs)
    print ("Counts matrix")
    print_matrix (motifs.counts)
    print ("PWM")
    print_matrix (motifs.pwm)
    print ("Sequence alphabet")
    print(motifs.alphabet)
    
    print(motifs.probability_sequence("AAACCT"))
    print(motifs.probability_sequence("ATACAG"))
    print(motifs.most_probable_sequence("CTATAAACCTTACATC"))
    
    for s in lseqs:
        print (s)
    print ("Consensus sequence")
    print(motifs.consensus())
    print ("Masked Consensus sequence")
    print(motifs.masked_consensus())

    s1 = MySeq("TAAAGTTATGA")
    s2 = MySeq("ATGACACGTG")
    s3 = MySeq("TTTGGGTAT")

    newmotif = motifs.create_motif([s1,s2,s3])
    print_matrix(newmotif.counts)
    print(newmotif.most_probable_sequence("AAAACT"))
    

if __name__ == '__main__':
    test()
