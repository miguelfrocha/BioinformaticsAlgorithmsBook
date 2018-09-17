# -*- coding: utf-8 -*-


def read_str(filename):
    f = open(filename)
    s = f.readline() ## assuming string is in first line
    return s

def convert_upper(s):
    return s.upper()

filename = input("Insert filename: ")
seq = read_str(filename)
seq1 = convert_upper(seq)
print("String in capital letters:", seq1)
