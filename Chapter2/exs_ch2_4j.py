# -*- coding: utf-8 -*-

def order_list_v1(l):
    lord = sorted(l, reverse = True)
    return lord

def order_list_v2(l):
    lord = []
    for x in l:
        i = 0
        while i < len(lord) and x < lord[i]:        
            i += 1
        lord.insert(i, x)
    return lord
          
# main
    
x = int(input("Value: "))
l = []
while x != 0:
    l.append(x)
    x = int(input("Value: "))
    
#lord = order_list_v1(l)
lord = order_list_v2(l)
print("Ordered list:", lord)