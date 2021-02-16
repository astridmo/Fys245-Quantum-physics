#!python
# -*- coding: utf-8 -*-

"""
...information about the code...
"""

__author__ = 'Astrid Moum'
__email__ = 'astridmo@nmbu.no'

import math as m
import cmath


sigma = 0.01
t = 4
h = 6.63 * 10**(-34)  # Plancks constant
mass = 0.0000000008
x = 2
v = 10
p0 = mass * v
comp_var = complex(0, (t*h)/mass)
k0 = 3
w0 = 4
#phi = (m.sqrt(sigma)/m.pi**(1/4) * m.sqrt(sigma**2 + complex(0, (t*h)/m))) #* m.exp(-(x-(p0*t))/m)
phi = (m.sqrt(sigma)/(m.pi**(1/4) * cmath.sqrt(sigma**2 + comp_var))) * cmath.exp((-(x-(p0*t)/mass)**2)/(2 * sigma**2 * comp_var)) * cmath.exp(complex(0,(k0*x - w0*t)))
print(phi)

t = 1
h = 2
m = 4
print(cmath.sqrt(sigma**2 + complex(0, (t*h)/m)))
#print(math.sqrt(3))