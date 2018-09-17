# -*- coding: utf-8 -*-

from MyGraph import MyGraph

class DeBruijnGraph (MyGraph):
    
    def __init__(self, frags):
        MyGraph.__init__(self, {})
        self.create_deBruijn_graph(frags)

    def add_edge(self, o, d):
        if o not in self.graph.keys():
            self.add_vertex(o)
        if d not in self.graph.keys():
            self.add_vertex(d)
        self.graph[o].append(d)

    def in_degree(self, v):
        res = 0
        for k in self.graph.keys(): 
            if v in self.graph[k]: 
                res += self.graph[k].count(v)
        return res

    def create_deBruijn_graph(self, frags):
        for seq in frags:
            suf = suffix(seq)
            self.add_vertex(suf)
            pref = prefix(seq)
            self.add_vertex(pref)
            self.add_edge(pref, suf)

    def seq_from_path(self, path):
        seq = path[0]
        for i in range(1,len(path)):
            nxt = path[i]
            seq += nxt[-1]
        return seq 
    
def suffix (seq): 
    return seq[1:]
    
def prefix(seq):
    return seq[:-1]

def composition(k, seq):
    res = []
    for i in range(len(seq)-k+1):
        res.append(seq[i:i+k])
    res.sort()
    return res

## exercise 2 - chapter 15
def tries_DeBruijn(frags):
    # try with original size
    dbgr = DeBruijnGraph(frags)
    nb = dbgr.check_nearly_balanced_graph()
    if (nb[0] is not None and nb[1] is not None):
        p = dbgr.eulerian_path()
        return dbgr.seq_from_path(p)
    
    k = len(frags[0])  # assuming all of the same size (not tested)
    while (k >= 2):
        nfrags = []
        for f in frags:
            nf = composition(k, f)
            nfrags.extend(nf)
        dbgr = DeBruijnGraph(nfrags)
        nb = dbgr.check_nearly_balanced_graph()
        if (nb[0] is not None and nb[1] is not None):
            p = dbgr.eulerian_path()
            return dbgr.seq_from_path(p)
        else: k -= 1
    return None


def test6():
    frags = [ "ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    dbgr = DeBruijnGraph(frags)
    dbgr.print_graph()
    print (dbgr.check_nearly_balanced_graph())
    print (dbgr.eulerian_path())


def test7():
    orig_sequence = "ATGCAATGGTCTG"
    frags = composition(3, orig_sequence)
    dbgr = DeBruijnGraph(frags)
    dbgr.print_graph()
    print (dbgr.check_nearly_balanced_graph())
    p= dbgr.eulerian_path()
    print (p)
    print (dbgr.seq_from_path(p))

def test8():
    frags = [ "ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA", "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    print (tries_DeBruijn(frags))

test6()
print()
test7()
print()
test8()
    