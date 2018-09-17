# -*- coding: utf-8 -*-
    
def sum_int_values_interval(lower, upper):
    ## assuming closed intervals
    s = 0
    for i in range(lower, upper + 1):
        s += i
    return s

a = int(input("Lower limit of interval: "))
b = int(input("Upper limit of interval: "))

s = sum_int_values_interval(a, b)

print("Sum of integer values in the interval", s)