# -*- coding: utf-8 -*-

x = int(input("Value: "))
lord = []

while x != 0:
    i = 0
    while i < len(lord) and x < lord[i]:        
        i += 1
    lord.insert(i, x)
    x = int(input("Value: "))
          
print("Ordered list:", lord)
