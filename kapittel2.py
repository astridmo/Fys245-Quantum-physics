#!python
# -*- coding: utf-8 -*-

"""
...information about the code...
"""

__author__ = 'Astrid Moum'
__email__ = 'astridmo@nmbu.no'

import math as m
import cmath
import numpy as np

h = 6.63 * 10 ** (-34)  # Plancks constant


def phi(sigma, t, x, mass, v):
    """Calculate wave function"""
    p0 = mass * v
    comp_var = complex(0, (t * h) / mass)


    phi = (m.sqrt(sigma) / (m.pi ** (1 / 4) * cmath.sqrt(sigma ** 2 + comp_var))) * cmath.exp(
        (-(x - (p0 * t) / mass) ** 2) / (2 * sigma ** 2 + comp_var)) * cmath.exp(complex(0, (k0 * x - w0 * t)))

    return phi

sigma = 0.01 * 10 ** (-2)
t = 5
mass = 9.1094 * 10 ** (-31)  # Mass of electron
x = 2
v = 10
k0 = 3
w0 = 4

print(phi(sigma, t, x, mass, v))
