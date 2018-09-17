# -*- coding: utf-8 -*-

seq=input("Insert sequence: ")

seq = seq.upper()

if seq == seq[::-1]:
    print("Palindrome")
else:
    print("Not palindrome")
    
## check exercise 4e for further alternatives