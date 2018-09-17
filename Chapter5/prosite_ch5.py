
def find_zync_finger(seq):
    from re import search
    regexp = "C.H.[LIVMFY]C.{2}C[LIVMYA]"
    mo = search(regexp, seq)
    if (mo != None):
        return mo.span()[0]
    else:
        return -1
    
def find_prosite(seq, profile):
    from re import search
    regexp = profile.replace("-","")
    regexp = regexp.replace("x",".")
    regexp = regexp.replace("(","{")
    regexp = regexp.replace(")","}")
    mo = search(regexp, seq)
    if (mo != None):
        return mo.span()[0]
    else:
        return -1
    
def test():
    seq = "HKMMLASCKHLLCLKCIVKLG"
    print(find_zync_finger(seq))
    print(find_prosite(seq,"C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"))
    
test()