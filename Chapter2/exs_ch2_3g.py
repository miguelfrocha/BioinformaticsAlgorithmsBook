# -*- coding: utf-8 -*-
    
a = int(input("Insert lower limit of first interval: "))
b = int(input("Insert upper limit of first interval: "))
c = int(input("Insert lower limit of second interval: "))
d = int(input("Insert upper limit of second interval: "))

inf_int = max(a, c)
sup_int = min(b, d)

inf_union = min(a,c)
sup_union = max(b,d)

if inf_int <= sup_int:
    inters = [inf_int, sup_int]
else: inters = []
print("Intersection of [{}, {}] and [{}, {}]: {}".format(a,b,c,d,inters))

if c <= b:
    union = [inf_union, sup_union]
    print("Union of [{}, {}] and [{}, {}]: {}".format(a,b,c,d,union))
else:
    print("Union of [{}, {}] and [{}, {}]: [{},{}] U [{},{}]".format(a,b,c,d,a,b,c,d))