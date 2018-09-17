# -*- coding: utf-8 -*-

a = float(input("Insert first value: "))
b = float(input("Insert second value: "))
c = float(input("Insert third value: "))

if a > b:
    if a > c: l = a
    else: l = c
    if b > c: s = c
    else: s = b
else:
    if b > c: l = b
    else: l = c
    if a > c: s = c
    else: s = a 

print("Largest value:", l)
print("Smallest value:", s)