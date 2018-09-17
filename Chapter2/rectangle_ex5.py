#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math

class Rectangle:
    
    def __init__(self, height, width):
        self.height = height
        self.width = width
        
    def __str__(self):
        return ("H:" + str(self.height) + " W:" + str(self.width))
        
    def area(self):
        return self.height * self.width
    
    def perimeter(self):
        return 2*self.height + 2*self.width
    
    def diagonal(self):
        return math.sqrt(self.height ** 2 + self.width ** 2)

class Square (Rectangle):
    
    def __init__(self, side):
        self.height = side
        self.width = side
        
r1 = Rectangle(3,4)
print(r1)
print("Area:", r1.area())
print("Perimeter:", r1.perimeter())
print("Diagnonal", r1.diagonal())

print()

s1 = Square(10)
print(s1)
print("Area:", s1.area())
print("Perimeter:", s1.perimeter())
print("Diagnonal", s1.diagonal())

