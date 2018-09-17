## general purpose
def find_pattern_re (seq, pat):
    from re import search
    mo = search(pat, seq)
    if (mo != None):
        return mo.span()[0]
    else:
        return -1

def find_all_occurrences_re (seq, pat):
    from re import finditer
    mos = finditer(pat, seq)
    res = []
    for x in mos:
        res.append(x.span()[0])
    return res

def find_all_overlap(seq, pat):
    return find_all_occurrences_re(seq, "(?="+pat+")")


def test():
    seq = input("Input sequence:")
    pat = input("Input pattern (as a regular expression):")
    
    res = find_pattern_re(seq, pat)
    if res >= 0:
        print("Pattern found in position: ", res)
    else:  print("Pattern not found")
    
    all_res = find_all_occurrences_re(seq, pat)
    if len(all_res) > 0:
        print("Pattern found in positions: ", all_res)
    else:  print("Pattern not found")
        
    all_ov = find_all_overlap(seq, pat)
    if len(all_ov) > 0:
        print("Pattern found in positions (overlap): ", all_ov)
    else: 
        print("Pattern not found")


test()