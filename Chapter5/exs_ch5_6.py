# -*- coding: utf-8 -*-

def handle_uniprot_header (header):
    from re import search

    dic = {}
    er = ">.*\|(.*)\|(\S*)(.*)OS=(.*)GN=(.*)PE=(.*)SV=(.*)"
    res =  search(er, header)
    er1 = ">.*\|(.*)\|(\S*)(.*)OS=(.*)PE=(.*)SV=(.*)"
    res1 = search(er1, header)
    if res != None:
        idprot = res.group(1)
        entry = res.group(2)
        desc = res.group(3)
        organism = res.group(4)
        gene = res.group(5)    
        pe = res.group(6)    
        sv = res.group(7)    
    elif res1 != None:
        idprot = res1.group(1)
        entry = res1.group(2)
        desc = res1.group(3)
        organism = res1.group(4)
        gene = None
        pe = res1.group(5)    
        sv = res1.group(6)        
    else: return None
    
    dic["id"] = idprot
    dic["Entry"] = entry
    dic["Protein"] = desc
    dic["OS"] = organism
    dic["PE"] = pe
    dic["SV"] = sv
    dic["GN"] = gene
    
    return dic
    
    
print ( handle_uniprot_header(">sp|P27747|ACOX_RALEH Acetoin catabolism protein X OS=Ralstonia eutropha (strain ATCC 17699 / H16 / DSM 428 / Stanier 337) GN=acoX PE=4 SV=2") )
print ( handle_uniprot_header(">sp|P27747|ACOX_RALEH Acetoin catabolism protein X OS=Ralstonia eutropha (strain ATCC 17699 / H16 / DSM 428 / Stanier 337) PE=4 SV=2") )