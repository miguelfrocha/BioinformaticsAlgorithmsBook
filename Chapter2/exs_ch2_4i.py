# -*- coding: utf-8 -*-

def sum_list(l):
    s = 0
    for e in l:
        s += e
    return s

def mean_list(l):
    s = sum_list(l)
    return s / len(l)

def max_list(l):
    m = l[0]
    for e in l:
        if e > m:
            m = e
    return m


l = []
x = int(input("Value: "))
l.append(x)
while x != 0:
    x = int(input("Value: "))
    l.append(x)
    
print("Sum: {}\tMean: {}\tMaximum: {}".format(sum_list(l), mean_list(l), max_list(l)))