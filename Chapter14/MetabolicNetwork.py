# -*- coding: utf-8 -*-
"""
Created on Thu May 18 20:01:45 2017

@author: miguelrocha
"""

from MyGraph import MyGraph

class MetabolicNetwork (MyGraph):
    
    def __init__(self, network_type = "metabolite-reaction", split_rev = False):
        MyGraph.__init__(self, {})
        self.net_type = network_type
        self.node_types = {}
        if network_type == "metabolite-reaction":
            self.node_types["metabolite"] = []
            self.node_types["reaction"] = []
        self.split_rev =  split_rev
    
    def add_vertex_type(self, v, nodetype):
        self.add_vertex(v)
        self.node_types[nodetype].append(v)
    
    def get_nodes_type(self, node_type):
        if node_type in self.node_types:
            return self.node_types[node_type]
        else: return None
    
    def load_from_file(self, filename):
        rf = open(filename)
        gmr = MetabolicNetwork("metabolite-reaction")
        for line in rf:
            if ":" in line:
                tokens = line.split(":")
                reac_id = tokens[0].strip()
                gmr.add_vertex_type(reac_id, "reaction")
                rline = tokens[1]
            else: raise Exception("Invalid line:")                
            if "<=>" in rline:
                left, right = rline.split("<=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    if self.split_rev:
                        gmr.add_vertex_type(reac_id+"_b", "reaction")
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id+"_b", met_id)
                    else:
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id, met_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    if self.split_rev:
                        gmr.add_edge(met_id, reac_id+"_b")
                        gmr.add_edge(reac_id, met_id)
                    else:
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id, met_id)
            elif "=>" in line:
                left, right = rline.split("=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    gmr.add_edge(met_id, reac_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    gmr.add_edge(reac_id, met_id)
            else: raise Exception("Invalid line:")    

        
        if self.net_type == "metabolite-reaction": 
            self.graph = gmr.graph
            self.node_types = gmr.node_types
        elif self.net_type == "metabolite-metabolite":
            self.convert_metabolite_net(gmr)
        elif self.net_type == "reaction-reaction": 
            self.convert_reaction_graph(gmr)
        else: self.graph = {}
        
        
    def convert_metabolite_net(self, gmr):
        for m in gmr.node_types["metabolite"]:
            self.add_vertex(m)
            sucs = gmr.get_successors(m)
            for s in sucs:
                sucs_r = gmr.get_successors(s)
                for s2 in sucs_r:
                    if m != s2: 
                        self.add_edge(m, s2)

        
    def convert_reaction_graph(self, gmr): 
        for r in gmr.node_types["reaction"]:
            self.add_vertex(r)
            sucs = gmr.get_successors(r)
            for s in sucs:
                sucs_r = gmr.get_successors(s)
                for s2 in sucs_r:
                    if r != s2: self.add_edge(r, s2)
      
    ## assessing metabolic potential      
      
    def active_reactions(self, active_metabolites):
        if self.net_type != "metabolite-reaction" or not self.split_rev:
            return None
        res = []
        for v in self.node_types['reaction']:
            preds = set(self.get_predecessors(v))
            if len(preds)>0 and preds.issubset(set(active_metabolites)):
                res.append(v)
        return res
        
    def produced_metabolites(self, active_reactions):
        res = []
        for r in active_reactions:
            sucs = self.get_successors(r)
            for s in sucs:
                if s not in res: res.append(s)
        return res
        
    def all_produced_metabolites(self, initial_metabolites):
        mets = initial_metabolites
        cont = True
        while cont:
            cont = False
            reacs = self.active_reactions(mets)
            new_mets = self.produced_metabolites(reacs)
            for nm in new_mets: 
                if nm not in mets: 
                    mets.append(nm)
                    cont = True
        return mets

    ## exercise 1 of chapter 14
    def final_metabolites(self):
        res = []
        for v in self.graph.keys():
            if v[0] == "M":
                if len(self.get_predecessors(v) ) > 0:
                    if self.get_successors(v) == []:
                        res.append(v)
        return res
    
    ## exercise 2 of chapter 14
    def shortest_path_product(self, initial_metabolites, target_product):
        if target_product in initial_metabolites: 
            return []
        metabs = {}
        for m in initial_metabolites: 
            metabs[m] = []
        reacs = self.active_reactions(initial_metabolites)
        cont = True
        while cont:
            cont = False
            for r in reacs:
                sucs = self.get_successors(r)
                preds = self.get_predecessors(r)
                for s in sucs:
                    if s not in metabs: 
                        previous = []
                        for p in preds:
                            for rr in metabs[p]:
                                if rr not in previous: previous.append(rr)
                        metabs[s] = previous + [r]
                        if s == target_product: 
                            return metabs[s]
                        cont = True
            if cont: 
                reacs = self.active_reactions(metabs.keys())
        return None
    
def test1():
    m = MetabolicNetwork("metabolite-reaction")
    m.add_vertex_type("R1","reaction")
    m.add_vertex_type("R2","reaction")
    m.add_vertex_type("R3","reaction")
    m.add_vertex_type("M1","metabolite")
    m.add_vertex_type("M2","metabolite")
    m.add_vertex_type("M3","metabolite")
    m.add_vertex_type("M4","metabolite")
    m.add_vertex_type("M5","metabolite")
    m.add_vertex_type("M6","metabolite")
    m.add_edge("M1","R1")
    m.add_edge("M2","R1")
    m.add_edge("R1","M3")
    m.add_edge("R1","M4")
    m.add_edge("M4","R2")
    m.add_edge("M6","R2")
    m.add_edge("R2","M3")
    m.add_edge("M4","R3")
    m.add_edge("M5","R3")
    m.add_edge("R3","M6")
    m.add_edge("R3","M4")
    m.add_edge("R3","M5")
    m.add_edge("M6","R3")
    m.print_graph()
    print("Reactions: ", m.get_nodes_type("reaction") )
    print("Metabolites: ", m.get_nodes_type("metabolite") )
    
    ## ex 1
    print(m.final_metabolites())

        
def test2():

    print("metabolite-reaction network:")
    mrn = MetabolicNetwork("metabolite-reaction")
    mrn.load_from_file("example-net.txt")
    mrn.print_graph()
    print("Reactions: ", mrn.get_nodes_type("reaction") )
    print("Metabolites: ", mrn.get_nodes_type("metabolite") )
    print()
    
    print("metabolite-metabolite network:")
    mmn = MetabolicNetwork("metabolite-metabolite")
    mmn.load_from_file("example-net.txt")
    mmn.print_graph()
    print()
    
    print("reaction-reaction network:")
    rrn = MetabolicNetwork("reaction-reaction")
    rrn.load_from_file("example-net.txt")
    rrn.print_graph()
    print()
    
    print("metabolite-reaction network (splitting reversible):")
    mrsn = MetabolicNetwork("metabolite-reaction", True)
    mrsn.load_from_file("example-net.txt")
    mrsn.print_graph()
    print()
    
    print("reaction-reaction network (splitting reversible):")
    rrsn = MetabolicNetwork("reaction-reaction", True)
    rrsn.load_from_file("example-net.txt")
    rrsn.print_graph()
    print()
    
    print(mmn.mean_degree("out")) 
    print(mmn.prob_degree("out"))
    
    print(mmn.mean_distances())
    print(mrn.mean_distances())
    
    print(mmn.all_clustering_coefs())
    print(mmn.mean_clustering_coef())
    print(mmn.mean_clustering_perdegree())
    
    print(mmn.highest_degrees(top = 3)) 
    print(mmn.highest_closeness(top = 3))
    
    print(mmn.betweenness_centrality("M5"))
  
def test3():
    print("metabolite-reaction network:")
    ec_mrn = MetabolicNetwork("metabolite-reaction")
    ec_mrn.load_from_file("ecoli.txt")
    print("Size:", ec_mrn.size())
    
    print("Mean degree: ", ec_mrn.mean_degree("inout")) 
    pd = ec_mrn.prob_degree("inout")
    pdo = sorted( list (pd.items()), key=lambda x : x[0]) 
    print("Histogram of degree probabilities")
    for (x,y) in pdo: print(x, "->", y)        
    print("Mean distance (M-R): ", ec_mrn.mean_distances())
    print()    
    
    print("metabolite-metabolite network:")
    ec_mmn = MetabolicNetwork("metabolite-metabolite")
    ec_mmn.load_from_file("ecoli.txt")
    print(ec_mmn.size())
    
    print("Mean degree: ", ec_mmn.mean_degree("inout")) 
    pd = ec_mmn.prob_degree("inout")
    pdo = sorted(list(pd.items()), key=lambda x : x[0])
    print("Histogram of degree probabilities")
    for (x,y) in pdo: print(x, "->", y)    
    
    print("Mean distance (M-M): ", ec_mmn.mean_distances()) 

    print(ec_mmn.mean_clustering_coef())
    cc = ec_mmn.mean_clustering_perdegree()
    cco = sorted(list(cc.items()), key=lambda x : x[0])
    print("Clustering coefficients per degree")
    for (x,y) in cco: print(x, "->", y)
        
    print(ec_mmn.highest_degrees(top = 20))
    print(ec_mmn.highest_closeness(top = 20))
    print()    
    
    print("reaction-reaction network:")
    ec_rrn = MetabolicNetwork("reaction-reaction")
    ec_rrn.load_from_file("ecoli.txt")
    print(ec_rrn.size())
    
    print("Mean degree: ", ec_rrn.mean_degree("inout")) 
    pd = ec_rrn.prob_degree("inout")
    pdo = sorted(list(pd.items()), key=lambda x : x[0])
    print("Histogram of degree probabilities")
    for (x,y) in pdo: print(x, "->", y)        
    
  
def test4():
    mrsn = MetabolicNetwork("metabolite-reaction", True)
    mrsn.load_from_file("example-net.txt")
    mrsn.print_graph()
    
    print(mrsn.produced_metabolites(["R1"]))
    print(mrsn.active_reactions(["M1","M2"]))
    print(mrsn.all_produced_metabolites(["M1","M2"]))
    print(mrsn.all_produced_metabolites(["M6"]))
    
    ## ex. 2
    print(mrsn.shortest_path_product(["M1","M2"], "M4"))
    print(mrsn.shortest_path_product(["M6"], "M3"))

test1()
print()
test2()
print()
#test3()
print()
test4()
