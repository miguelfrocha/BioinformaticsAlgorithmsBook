from Bio import Phylo

tree = Phylo.read("simple.dnd", "newick")
print(tree)

Phylo.draw_ascii(tree)

tree2 = Phylo.read("int_node_labels.nwk", "newick") 
Phylo.draw_ascii(tree2)

Phylo.convert("int_node_labels.nwk", "newick", "tree.xml", "phyloxml")
trees = Phylo.parse("tree.xml", "phyloxml") 
for t in trees: print(t)
    
from Bio.Phylo.PhyloXML import Phylogeny

treep = Phylogeny.from_tree(tree)
Phylo.draw(treep)

treep.root.color = "gray"
mrca = treep.common_ancestor({"name": "E"}, {"name": "F"})
mrca.color = "salmon"
treep.clade[0, 1].color = "blue"
Phylo.draw(treep)
 
