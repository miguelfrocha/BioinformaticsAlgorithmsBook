# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 02:22:38 2017

@author: miguelrocha
"""

from Bio.Seq import Seq 
my_seq = Seq("ATAGAGAAATCGCTGC")

print(my_seq) 
print(my_seq.alphabet)

from Bio.Alphabet import IUPAC 
my_seq = Seq("ATAGAGAAATCGCTGC", IUPAC.unambiguous_dna) 
print(my_seq)

my_prot = Seq("MJKLKVERSVVMSVLP", IUPAC.protein)
print(my_prot)

print(IUPAC.unambiguous_dna.letters)
print(IUPAC.ambiguous_dna.letters)
print(IUPAC.IUPACProtein.letters)
print(IUPAC.ExtendedIUPACProtein.letters)

for i in my_seq: print(i)
print(len(my_seq))
print(my_seq[2:4])
print(my_seq.count("G"))
print("GAGA" in my_seq)
print(my_seq.find("ATC"))

seq1 = Seq("MEVRNAKSLV", IUPAC.protein)
seq2 = Seq("GHERWKY", IUPAC.protein)
print(seq1+seq2)

from Bio.Alphabet import generic_nucleotide
nuc_seq = Seq("ATAGAGAAATCGCTGC", generic_nucleotide) 
dna_seq = Seq("TGATAGAACGT", IUPAC.unambiguous_dna) 
print(nuc_seq + dna_seq )

coding_dna = Seq("ATGAAGGCCATTGTAATGGGCCGC", IUPAC.unambiguous_dna) 
template_dna = coding_dna.reverse_complement() 
print( template_dna )
messenger_rna = coding_dna.transcribe() 
print(messenger_rna)

rna_seq = Seq('AUGCGUUUAACU', IUPAC.unambiguous_rna)
print(rna_seq.translate())
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
print(coding_dna.translate())
print(coding_dna.translate(table="Vertebrate Mitochondrial") )

from Bio.Data import CodonTable 
standard_table = CodonTable.unambiguous_dna_by_name["Standard"] 
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
print (standard_table)
print(mito_table.stop_codons) 
['TAA', 'TAG', 'AGA', 'AGG'] 
print(mito_table.start_codons) 
['ATT', 'ATC', 'ATA', 'ATG', 'GTG'] 
print(mito_table.forward_table["ATA"])
print(standard_table.forward_table["ATA"]) 


## Sect. sequence annotations

from Bio import SeqFeature
start_pos = SeqFeature.AfterPosition(5)
end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
my_location = SeqFeature.FeatureLocation(start_pos, end_pos)
print(my_location)
print(int(my_location.start) )
print(int(my_location.end) )

example_parent = Seq("ACCGAGACGGCAAAGGCTAGCATAGGTATGAGACTT")
from Bio.SeqFeature import SeqFeature, FeatureLocation
example_feature = SeqFeature(FeatureLocation(5, 18), type="gene", strand=-1)
feature_seq = example_feature.extract(example_parent)
print(feature_seq)

from Bio.Seq import Seq
seq = Seq("ATGAATGATAGCTGAT")
from Bio.SeqRecord import SeqRecord
seq_rec = SeqRecord(seq)
seq_rec.id = "ABC12345"
seq_rec.description = "My own sequence."
seq_rec.annotations["role"] = "unknown"
print(seq_rec)



from Bio import SeqIO 
record = SeqIO.read("NC_005816.fna", "fasta") 
print(record)
print(len(record.seq))
print(record.id)
print(record.description)
print(record.dbxrefs)
print(record.annotations)
print(record.features)

print("###### Formato Genbank")

record = SeqIO.read("NC_005816.gb", "genbank")
print(record.seq)
print(record.id)
print(record.name)
print(record.description)
print(record.dbxrefs)
print(len(record.annotations) )
print(record.annotations)
print(record.annotations["source"] )
print(record.annotations["taxonomy"] )
print(record.annotations["date"] )
print(record.annotations["gi"] )

print(record.letter_annotations)
print(record.features)

feat_genes = []
for i in range(len(record.features)):
    if record.features[i].type == "gene":
        feat_genes.append(record.features[i])
print(len(feat_genes))
for f in feat_genes: print(f.qualifiers['locus_tag'], f.strand, f.location)

for f in feat_genes: print(f.extract(record.seq).translate(table="Bacterial", cds=True))
    
feat_cds = []
for i in range(len(record.features)):
    if record.features[i].type == "CDS":
        feat_cds.append(record.features[i])    
print(len(feat_cds))

for (f1,f2) in zip(feat_genes, feat_cds):
    translated = f1.extract(record.seq).translate(table="Bacterial", cds=True)
    cdsprot = f2.qualifiers['translation']
    print(translated == cdsprot[0])

all_species = []
for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"): 
    print (seq_record.description)
    all_species.append(seq_record.annotations["organism"]) 
print(all_species)

from Bio import Entrez 
from Bio import SeqIO 
Entrez.email = "example@gmail.com" 
handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="6273291, 6273290, 6273289") 
for seq_record in SeqIO.parse(handle, "gb"):
	print (seq_record.id, seq_record.description[:100], "...")
	print ("Sequence length: ", len(seq_record)) 
handle.close() 