# -*- coding: utf-8 -*-

s = 0
n = 0.0

x = int(input("Value: "))
m = x
while x != 0:
    s += x
    n += 1
    if m < x:
        m = x
    x = int(input("Value: "))
        
print("Sum: {}\tMean: {}\tMaximum: {}".format(s, s / n, m))