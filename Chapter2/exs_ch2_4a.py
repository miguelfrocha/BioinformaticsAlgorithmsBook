#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def fahrenheit(celsius):
    f = 32 + 9 * celsius / 5.0
    return f

def test_fahrenheit():
    t = int(input("Temperature - Celsius: "))
    f = fahrenheit(t)
    print("Temperature in Fahrenheit:", f)
    

test_fahrenheit()
