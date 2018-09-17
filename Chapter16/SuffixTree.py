# -*- coding: utf-8 -*-

class SuffixTree:
    
    def __init__(self):
        self.nodes = { 0:(-1,{}) } # root node
        self.num = 0
    
    def print_tree(self):
        for k in self.nodes.keys():
            if self.nodes[k][0] < 0:
                print (k, "->", self.nodes[k][1]) 
            else:
                print (k, ":", self.nodes[k][0])
                
    def add_node(self, origin, symbol, leafnum = -1):
        self.num += 1
        self.nodes[origin][1][symbol] = self.num
        self.nodes[self.num] = (leafnum,{})
        
    def add_suffix(self, p, sufnum):
        pos = 0
        node = 0
        while pos < len(p):
            if p[pos] not in self.nodes[node][1].keys():
                if pos == len(p)-1:
                    self.add_node(node, p[pos], sufnum)
                else:
                    self.add_node(node, p[pos])     
            node = self.nodes[node][1][p[pos]]
            pos += 1
    
    def suffix_tree_from_seq(self, text):
        t = text+"$"
        for i in range(len(t)):
            self.add_suffix(t[i:], i)
            
    def find_pattern(self, pattern):
        pos = 0
        node = 0
        for pos in range(len(pattern)):
            if pattern[pos] in self.nodes[node][1].keys():
                node = self.nodes[node][1][pattern[pos]]
                pos += 1
            else: return None
        return self.get_leafes_below(node)
        

    def get_leafes_below(self, node):
        res = []
        if self.nodes[node][0] >=0: 
            res.append(self.nodes[node][0])            
        else:
            for k in self.nodes[node][1].keys():
                newnode = self.nodes[node][1][k]
                leafes = self.get_leafes_below(newnode)
                res.extend(leafes)
        return res
     
    ## exercise 1- chapter 16
    def find_parent(self, node):
        for k in self.nodes.keys():  
            if node in self.nodes[k][1].values():
                return k
        return None    
    
    def path_to_node (self, node):
        current = node
        res = ""
        while current != 0: ## root
            parent = self.find_parent(current)
            keys_parent = list(self.nodes[parent][1].keys())
            for kp in keys_parent:
                if self.nodes[parent][1][kp] == current: 
                    symb_par = kp
                    break
            res = symb_par + res
            current = parent
        return res
    
    def find_nodes_sizeK(self, k):
        active_nodes = [0]
        while (k>0):
            next_nodes = []
            for n in active_nodes:
                next_nodes.extend(list(self.nodes[n][1].values()))
            k = k - 1
            active_nodes = next_nodes
        return active_nodes
    
    def repeats(self, k, ocs):
        res = []
        active_nodes = self.find_nodes_sizeK(k)
        for n in active_nodes:
            ocs_pat = self.get_leafes_below(n)
            if len(ocs_pat) >= ocs:
                res.append(self.path_to_node(n))
        return res

def test():
    seq = "TACTA"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq)
    st.print_tree()
    print (st.find_pattern("TA"))
    print (st.find_pattern("ACG"))

def test2():
    seq = "TACTA"
    st = SuffixTree()
    st.suffix_tree_from_seq(seq)
    print (st.find_pattern("TA"))
    print(st.repeats(2,2))

test()
print()
test2()
        
            
    
    