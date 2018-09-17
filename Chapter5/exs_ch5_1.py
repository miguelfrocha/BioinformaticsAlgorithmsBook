#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def repeats(seq, k):
    dic = {}
    for i in range(len(seq)-k+1):
        cod = seq[i:i+k]
        if cod in dic: 
            dic[cod] += 1 
        else: dic[cod] = 1
    res = {}
    for k in dic.keys():
        if dic[k] > 1: res[k] = dic[k] 
    return res

print(repeats("ATAGAGATAGGAAGA",4))
print(repeats("ATAGAGATAGGAAGA",3))