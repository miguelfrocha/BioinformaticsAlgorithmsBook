# -*- coding: utf-8 -*-

def first_prot(seqAA):
    inProt = False
    endProt = False
    prot = ""
    pos = 0
    while not endProt and pos < len(seqAA):
        if inProt:
            if seqAA[pos]=="_": endProt = True
            else: prot += seqAA[pos]
        else:
            if seqAA[pos]=="M":
                inProt=True
                prot = "M"
        pos += 1
    if endProt: 
        return len(prot)
    else: 
        return -1

## test the function
seqAA= input("Aminoacid sequence:")
lprot = first_prot(seqAA)
if lprot>0: print ("Protein size: ", lprot)
else: print("Protein not found")
