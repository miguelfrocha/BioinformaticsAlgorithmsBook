from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio import SeqIO

##### comment this part after running once
record = SeqIO.read(open("interl10.fasta"), format="fasta") 
print (len(record.seq))
result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"))
save_file = open("interl-blast.xml", "w")
save_file.write(result_handle.read()) 
save_file.close() 
result_handle.close()
#####

result_handle = open("interl-blast.xml")
blast_record = NCBIXML.read(result_handle)

print ("PARAMETERS:")
print ("Database: " + blast_record.database)
print ("Matrix: " + blast_record.matrix)
print ("Gap penalties: ", blast_record.gap_penalties)
   
print ("number hits: ", len(blast_record.alignments) )
first_alignment = blast_record.alignments[0]

print ("FIRST ALIGNMENT:")
print ("Accession: " + first_alignment.accession)
print ("Hit id: " + first_alignment.hit_id)
print ("Definition: " + first_alignment.hit_def)
print ("Alignment length: ", first_alignment.length)
print ("Number of HSPs: ", len(first_alignment.hsps))


hsp = first_alignment.hsps[0]
print ("E-value: ", hsp.expect)
print ("Score: ", hsp.score)
print ("Length: ", hsp.align_length)
print ("Identities: ", hsp.identities)
print ("Alignment of the HSP:")
print (hsp.query)
print (hsp.match)
print (hsp.sbjct)

print ("Top 10 alignments:")
for i in range(10):
    alignment = blast_record.alignments[i]
    print ("Accession: " + alignment.accession)
    print ("Definition: " + alignment.hit_def)
    for hsp in alignment.hsps:
        print ("E-value: ", hsp.expect)
    print()

import re
specs = []
for i in range(20):
    alignment = blast_record.alignments[i]
    definition = alignment.hit_def
    x = re.search("\[(.*?)\]", definition).group(1)
    specs.append(x)

print ("Organisms:")
for s in specs: print(s)
