
def count_bases (seq):
    dic = {}
    seqC = seq.upper()
    errors = 0
    for b in seqC:
        if b in "ACGT":
            if b in dic: dic[b] += 1
            else: dic[b] = 1
        else: errors += 1
    return dic, errors
    
def print_perc_dic (dic):
    sum_values = sum(dic.values())
    for k in sorted(dic.keys()):
        print(" %s ->" % k, " %3.2f" % (dic[k]*100.0/sum_values), "%")

## main program
seq = input("Input DNA sequence: ")
freqs, errors = count_bases(seq)
if errors > 0: 
    print ("Sequence is invalid with ", errors , "invalid characters")
else: print("Sequence is valid")
print("Frequencies of the valid characters:")
print_perc_dic (freqs)

