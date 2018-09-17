# -*- coding: utf-8 -*-


def palindrome_v1(seq):
    if seq == seq[::-1]:
        return True
    else:
        return False

def palindrome_v2(seq):
    pal = True
    i = 0
    j = len(seq) - 1
    while pal and i < j:
        if seq[i] != seq[j]:
            pal = False
        i += 1
        j -= 1
    return pal

seq=input("Insert sequence: ")
if palindrome_v2(seq):
    print("Palindrome")
else:
    print("Not palindrome")
