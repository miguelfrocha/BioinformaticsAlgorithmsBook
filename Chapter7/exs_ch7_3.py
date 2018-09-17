
from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio import SeqIO

seqrecord = SeqIO.read(open("O14727.fasta"), format="fasta") 
print (len(seqrecord.seq))
  
result_handle = NCBIWWW.qblast("blastp", "swissprot", seqrecord.format("fasta"))

record = NCBIXML.read(result_handle)

print ("PARAMETERS:")
print ("Database: " + record.database)
print ("Matrix: " + record.matrix)
print ("Gap penalties: ", record.gap_penalties)

nhits = len(record.alignments) 
print ("number hits: ", nhits)

res = []
for alignment in record.alignments:
    evalue = alignment.hsps[0].expect
    accession = alignment.accession
    leng = alignment.hsps[0].align_length
    res.append(accession + " - " + str(evalue) + " length:" + str(leng) )

print("E-values and length of alignments:")
for s in res: print(s)

result_handle2 = NCBIWWW.qblast("blastp", "swissprot", seqrecord.format("fasta"), entrez_query = "Saccharomyces cerevisiae[organism]" )

record2 = NCBIXML.read(result_handle2)

first_alignment = record2.alignments[0]

print ("Accession: " + first_alignment.accession)
print ("Hit id: " + first_alignment.hit_id)
print ("Definition: " + first_alignment.hit_def)

print ("Number of HSPs: ", len(first_alignment.hsps))

for hsp in first_alignment.hsps:
    print ("E-value: ", hsp.expect)
    print ("Length: ", hsp.align_length)
    print ("Identities: ", hsp.identities)
    print ("Query start: ", hsp.query_start)
    print ("Sbjct start: ", hsp.sbjct_start)
    print (hsp.query[0:90])
    print (hsp.match[0:90])
    print (hsp.sbjct[0:90])
    print ("")
    