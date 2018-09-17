# -*- coding: utf-8 -*-

class Automata:
    
    def __init__(self, alphabet, pattern):
        self.numstates = len(pattern) + 1
        self.alphabet = alphabet
        self.transition_table = {}
        self.build_transition_table(pattern)        
    
    def build_transition_table(self, pattern):
        for q in range(self.numstates):
            for a in self.alphabet:
                prefix = pattern[0:q] + a
                self.transition_table[(q,a)] = overlap(prefix, pattern)
       
    def print_automata(self):
        print ("States: " , self.numstates)
        print ("Alphabet: " , self.alphabet)
        print ("Transition table:")
        for k in self.transition_table.keys():
            print (k[0], ",", k[1], " -> ", self.transition_table[k])
         
    def next_state(self, current, symbol):
        return self.transition_table.get((current, symbol))
        
    def apply_seq(self, seq):
        q = 0
        res = [q]
        for c in seq:
            q = self.next_state(q, c)
            res.append(q)
        return res
        
    def occurences_pattern(self, text):
        q = 0 
        res = []
        for i in range(len(text)):
            q = self.next_state(q, text[i])
            if q == self.numstates-1: 
                res.append(i - self.numstates + 2)
        return res

def overlap(s1, s2):
    maxov = min(len(s1), len(s2))
    for i in range(maxov,0,-1):
        if s1[-i:] == s2[:i]: return i
    return 0
               
def test():
    auto = Automata("ACGT", "ACA")
    auto.print_automata()
    print (auto.apply_seq("CACATGACATG"))
    print (auto.occurences_pattern("CACATGACATG"))

test()
        

