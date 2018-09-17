# -*- coding: utf-8 -*-
    
a = int(input("Lower limit of interval: "))
b = int(input("Upper limit of interval: "))

## assuming closed intervals

s = 0
for i in range(a, b + 1):
    s += i

print("Sum of integer values in the interval", s)