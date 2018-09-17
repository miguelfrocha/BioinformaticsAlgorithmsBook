#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dotplots import dotplot

def largest_diagonal(mat):
    largest_size = 0
    begin_largest = (-1,-1)

    # testar diagonais acima da principal inc principal
    for col in range(0,len(mat[0])):
        cur_size = 0
        for k in range(len(mat[0])-col):
            if mat[k][col+k] == 1: 
                if cur_size == 0: ini = (k, col + k)
                cur_size += 1
            if mat[k][col+k] == 0 or k == len(mat[0])-col-1:
                if cur_size > 0: 
                    if cur_size > largest_size:
                        largest_size = cur_size
                        begin_largest = ini
                    cur_size = 0
    # testar diagonais abaixo da principal
    for lin in range(1,len(mat)):
        cur_size = 0
        for k in range(len(mat)-lin):
            if mat[lin+k][k] == 1: 
                if cur_size == 0: ini = (lin+k, k)
                cur_size += 1
            if mat[lin+k][k] == 0 or k == len(mat)-lin-1:
                if cur_size > 0: 
                    if cur_size > largest_size:
                        largest_size = cur_size
                        begin_largest = ini
                    cur_size = 0                          
    return largest_size, begin_largest[0], begin_largest[1]


def test_diag():
    s1 = "ATAATA"
    s2 = "ATAATA"
    mat1 = dotplot(s1, s2)
    print(largest_diagonal(mat1))    
    
    s1 = "CGATATAG"
    s2 = "TATATATT"
    mat2 = dotplot(s1, s2)
    print(largest_diagonal(mat2))
    
test_diag()