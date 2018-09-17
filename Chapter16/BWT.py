# -*- coding: utf-8 -*-

class BWT:
    
    def __init__(self, seq = "", buildsufarray = False):
        self.bwt = self.build_bwt(seq, buildsufarray) 
    
    def set_bwt(self, bw):
        self.bwt = bw

    def build_bwt(self, text, buildsufarray = False):
        ls = []
        for i in range(len(text)):
            ls.append(text[i:]+text[:i])
        ls.sort()
        res = ""
        for i in range(len(text)):
            res += ls[i][len(text)-1]
        if buildsufarray:
            self.sa = []
            for i in range(len(ls)):
                stpos = ls[i].index("$")
                self.sa.append(len(text)-stpos-1)
        print(self.sa)
        return res    
    
    def inverse_bwt(self):
        firstcol = self.get_first_col()
        res = ""
        c = "$" 
        occ = 1
        for i in range(len(self.bwt)):
            pos = find_ith_occ(self.bwt, c, occ)
            c = firstcol[pos]
            occ = 1
            k = pos-1
            while firstcol[k] == c and k >= 0: 
                occ += 1
                k -= 1
            res += c         
        return res        
 
    def get_first_col (self):
        firstcol = []
        for c in self.bwt: firstcol.append(c)
        firstcol.sort()  
        return firstcol
        
    def last_to_first(self):
        res = []
        firstcol = self.get_first_col()
        for i in range(len(firstcol)):
            c = self.bwt[i]
            ocs = self.bwt[:i].count(c) + 1
            res.append(find_ith_occ(firstcol, c, ocs))
        return res

    def bw_matching(self, patt):
        lf = self.last_to_first()
        res = []
        top = 0
        bottom = len(self.bwt)-1
        flag = True
        while flag and top <= bottom:
            if patt != "":
                symbol = patt[-1]
                patt = patt[:-1]
                lmat = self.bwt[top:(bottom+1)]
                if symbol in lmat:
                    topIndex = lmat.index(symbol) + top
                    bottomIndex = bottom - lmat[::-1].index(symbol)
                    top = lf[topIndex]
                    bottom = lf[bottomIndex]
                else: flag = False
            else: 
                for i in range(top, bottom+1): res.append(i)
                flag = False            
        return res        
 
    def bw_matching_pos(self, patt):
        res = []
        matches = self.bw_matching(patt)
        for m in matches:
            res.append(self.sa[m])
        res.sort()
        return res
 
# auxiliary
 
def find_ith_occ(l, elem, index):
    j,k = 0,0
    while k < index and j < len(l):
        if l[j] == elem:
            k = k +1
            if k == index: return j
        j += 1 
    return -1 


      
def test():
    seq = "TAGACAGAGA$"
    bw = BWT(seq, True)
    print (bw.bwt)
    print (bw.last_to_first())
    print (bw.bw_matching("AGA"))
    print("Suffix array:", bw.sa)
    print(bw.bw_matching_pos("AGA"))

def test2():
    bw = BWT()
    bw.set_bwt("ACG$GTAAAAC")
    print (bw.inverse_bwt())


test()