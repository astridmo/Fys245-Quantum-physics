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
    k0 = (2 * m.pi * p0) / h
    w0 = (2 * m.pi) / t
    phi = []

    for xi in x:
        phi.append(m.sqrt(sigma) / (m.pi ** (1 / 4) * cmath.sqrt(sigma ** 2 + (complex(0, h * t) / mass))) * cmath.exp(
            -0.5 * (xi - (p0 * t / mass)) ** 2 / (sigma ** 2 + (complex(0, h * t) / mass))) * cmath.exp(
            complex(0, k0 * xi - w0 * t)))

    return phi


def phi_squared(sigma, t, x, mass, v):
    """Calculate squared of wave function."""
    return np.square(phi(sigma, t, x, mass, v))


sigma = 0.01 * 10 ** (-2)
t = 9 * 10 ** (-5)
mass = 9.1094 * 10 ** (-31)  # Mass of electron
x = np.linspace(0, 10, 1000)
v = 10 ** 6 # m/s
k0 = 3
w0 = 4

# for xi in x:
#     print(phi(sigma, t, xi, mass, v))

# print(np.square(phi(sigma, t, x, mass, v)))
print(phi_squared(sigma, t, x, mass, v))
