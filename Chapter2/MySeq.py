# -*- coding: utf-8 -*-

class MySeq: 
    """ Class for biological sequences. """
    
    def __init__ (self, seq, seq_type = "DNA"): 
        self.seq = seq.upper()
        self.seq_type = seq_type

    def __len__(self):
        return len(self.seq)
    
    def __getitem__(self, n):
        return self.seq[n]

    def __getslice__(self, i, j):
        return self.seq[i:j]

    def __str__(self):
        return self.seq
    
    def print_sequence(self):
        print ("Sequence: " + self.seq)
    
    def get_seq_biotype (self):
        return self.seq_type
        
    def set_seq_biotype (self, bt):
        biotype = bt.upper()
        if biotype == "DNA" or biotype == "RNA" or biotype == "PROTEIN":
            self.seq_type = biotype
        else:
            print("Non biological sequence type!")    

    def show_info_seq (self):
        print ("Sequence: " + self.seq + " biotype: " + self.seq_type)
        
    def count_occurrences(self , seq_search): 
        return self.seq.count(seq_search)


s1 = MySeq("ATAATGATAGATAGATGAT")
## access attribute values
print (s1.seq)
print (s1.seq_type)
# calling methods
s1.print_sequence()
print (s1.get_seq_biotype())

s1.set_seq_biotype("time series")
s1.set_seq_biotype("dna")
print (s1.get_seq_biotype())

s1 = MySeq("MKKVSJEMSSVPYW", "PROTEIN")
print(s1)
print(len(s1))
print(s1[4])
print(s1[2:5])

