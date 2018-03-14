#!/usr/bin/env python3

h = 4.135e-15 # in eV
c = 3e8 # in m/s

def e_to_nu(e):
    return e/h

def e_to_lambda(e):
    return c / e_to_nu(e)


print(e_to_lambda(10.2 * 3 / 4))
print(e_to_lambda(10.2 * 3 / 4))
print(e_to_nu(10.2 * 1 / 2)/1e14)
