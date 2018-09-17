# -*- coding: utf-8 -*-

import math

def hypotenuse(s1, s2):
    h = math.sqrt(s1 ** 2 + s2 ** 2)
    return h

a = int(input("Length of side A: "))
b = int(input("Length of side B: "))

h = hypotenuse(a, b)

print("Length of hypotenuse:", h)
