# -*- coding: utf-8 -*-
from Bio.Blast import NCBIWWW 
from Bio import SeqIO
from Bio.Blast import NCBIXML 

##### comment this part after running first time 
## running blast
record = SeqIO.read(open("example_blast.fasta"), format="fasta") 
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))
## saving results
save_file = open("my_blast.xml", "w")
save_file.write(result_handle.read()) 
save_file.close() 
result_handle.close()
######

# analysing results

result_handle = open("my_blast.xml")

e_value_threshold = 0.001
blast_records = NCBIXML.parse(result_handle)
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < e_value_threshold:
                print('****Alignment****')
                print('sequence: ', alignment.title)
                print('length: ', alignment.length)
                print('e value: ', hsp.expect)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')