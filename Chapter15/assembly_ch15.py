# -*- coding: utf-8 -*-

from MyGraph import MyGraph

class OverlapGraph(MyGraph):
    
    def __init__(self, frags, reps = True):
        MyGraph.__init__(self, {})
        if reps: self.create_overlap_graph_with_reps(frags)
        else: self.create_overlap_graph(frags)
        self.reps = reps
        
    
    ## create overlap graph from list of sequences (fragments)
    def create_overlap_graph(self, frags):
        for seq in frags:
            self.add_vertex(seq)
        for seq in frags:
            suf = suffix(seq)
            for seq2 in frags:
                if prefix(seq2) == suf:
                    self.add_edge(seq, seq2)
    
    def create_overlap_graph_with_reps(self, frags):
        idnum = 1
        for seq in frags:
            self.add_vertex(seq+ "-" + str(idnum))
            idnum = idnum + 1
        idnum = 1
        for seq in frags:
            suf = suffix(seq)
            for seq2 in frags:
                if prefix(seq2) == suf:
                    for x in self.get_instances(seq2):
                        self.add_edge(seq+ "-" + str(idnum), x)    
            idnum = idnum + 1
    
    def get_instances(self, seq):
        res = []
        for k in self.graph.keys():
            if seq in k: res.append(k)
        return res
    
    def get_seq(self, node):
        if node not in self.graph.keys(): return None
        if self.reps: return node.split("-")[0]
        else: return node
    
    def seq_from_path(self, path):
        if not self.check_if_hamiltonian_path(path): return None
        seq = self.get_seq(path[0])
        for i in range(1,len(path)):
            nxt = self.get_seq(path[i])
            seq += nxt[-1]
        return seq    
   
                    
# auxiliary
def composition(k, seq):
    res = []
    for i in range(len(seq)-k+1):
        res.append(seq[i:i+k])
    res.sort()
    return res
    
def suffix (seq): 
    return seq[1:]
    
def prefix(seq):
    return seq[:-1]

  
# testing / mains
def test1():
    seq = "CAATCATGATG"
    k = 3
    print (composition(k, seq))
   
def test2():
    frags = ["ACC", "ATA", "CAT", "CCA", "TAA"]
    ovgr = OverlapGraph(frags, False)
    ovgr.print_graph()

def test3():
     frags = [ "ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
     ovgr = OverlapGraph(frags, True)
     ovgr.print_graph()
     path = ['ACC-2', 'CCA-8', 'CAT-5', 'ATG-3', 'TGG-13', 'GGC-10', 'GCA-9', 'CAT-6', 'ATT-4', 'TTT-15', 'TTC-14', 'TCA-12', 'CAT-7', 'ATA-1', 'TAA-11']
     print (ovgr.seq_from_path(path))
     print(ovgr.check_if_valid_path(path))
     print (ovgr.check_if_hamiltonian_path(path))
     

def test4():
    frags = [ "ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    ovgr = OverlapGraph(frags, True)
    ovgr.print_graph()
    path = ovgr.search_hamiltonian_path()
    print(path)
    print (ovgr.check_if_hamiltonian_path(path))
    print (ovgr.seq_from_path(path))

def test5():
    orig_sequence = "CAATCATGATGATGATC"
    frags = composition(3, orig_sequence)
    print (frags)
    ovgr = OverlapGraph(frags, True)
    ovgr.print_graph()
    path = ovgr.search_hamiltonian_path()
    print (path)
    print (ovgr.seq_from_path(path))
   
test1()
print()
test2()
print()
test3()
print()
test4()
print()
test5()
