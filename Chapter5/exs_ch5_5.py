# -*- coding: utf-8 -*-

def read_Fasta (filename):
    from re import sub, search

    res = []
    sequence = None
    info = None
    
    fh = open(filename)

    for line in fh:
        if search(">.*", line):
                if sequence is not None and info is not None and sequence != "":
                    res.append(sequence)
                info = line              
                sequence = ""
        else:
            if sequence is None: return None
            else: sequence += sub("\s","",line)

    if sequence is not None and info is not None and sequence != "":
                    res.append(sequence)
    fh.close()

    return res
    
print ( read_Fasta("exuniprot.fasta") )
print ( read_Fasta("NC_005816.fna") )