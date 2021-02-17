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
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy


h = 6.63 * 10 ** (-34)  # Plancks constant
sigma = 0.01 * 10 ** (-2)


def phi(sigma, t, x, mass, v):
    """Calculate wave function"""
    p0 = mass * v
    k0 = (2 * m.pi * p0) / h
    w0 = (m.pi * p0 ** 2) / (h * mass)
    phi = []

    for xi in x:
        phi.append(m.sqrt(sigma) / (m.pi ** (1 / 4) * cmath.sqrt(sigma ** 2 + (complex(0, h * t) / mass))) * cmath.exp(
            -0.5 * (xi - (p0 * t / mass)) ** 2 / (sigma ** 2 + (complex(0, h * t) / mass))) * cmath.exp(
            complex(0, k0 * xi - w0 * t)))

    return phi


def phi_squared(x, sigma, t, mass, v):
    """Calculate squared of wave function by using phi-function."""
    return np.square(phi(sigma, t, x, mass, v))


# def phi_squared2(x, sigma, t, mass, v):
#     """
#     Calculate squared of wave function.
#     NB! Tried to calculate by hand. Did not get it right (?)
#     """
#     p0 = mass * v
#     phi_squared = []
#     for xi in x:
#         phi_squared.append((-(xi-(p0*t)/mass)**2 * sigma**2) / (sigma**4 + ((t*h)/mass)**2))
#
#     return phi_squared

# def phi_squared_integral(x, sigma, t, mass, v):
#     """Calculation of pfi squared for use with integration"""
# https://stackoverflow.com/questions/5965583/use-scipy-integrate-quad-to-integrate-complex-numbers
#     p0 = mass * v
#     k0 = (2 * m.pi * p0) / h
#     w0 = (m.pi * p0 ** 2) / (h * mass)
#
#     def real_func(x):
#         return scipy.real(phi_squared_integral(x))
#
#     def imag_func(x):
#         return scipy.imag(phi_squared_integral(x))
#
#     return (m.sqrt(sigma) / (m.pi ** (1 / 4) * cmath.sqrt(sigma ** 2 + (complex(0, h * t) / mass))) * cmath.exp(
#             -0.5 * (x - (p0 * t / mass)) ** 2 / (sigma ** 2 + (complex(0, h * t) / mass))) * cmath.exp(
#             complex(0, k0 * x - w0 * t))) ** 2

# ==================
# Task a
# ==================


t1 = 0
t2 = 9 * 10 ** (-5)
m_elec = 9.1094 * 10 ** (-31)  # Mass of electron
x = np.linspace(-0.0003, 0.0003, 200)
v = 10 ** 6 # m/s
#x = 0.001
wave1 = phi_squared(x, sigma,t1, m_elec, v)
wave2 = phi_squared(x, sigma, t2, m_elec, v)
print(wave1)
print(wave2)
#
# wave = phi_squared(x, sigma, t1, m_elec, v)
#
plt.plot(x, wave1, label='electron t1')
plt.legend()
plt.show()
plt.plot(x, wave2, label='electron t2')
plt.legend()
plt.show()


# ====================
# Task b
# ====================

# Variables
m_helium = 6.64647 * 10 ** (-27)
v1 = 2.2 * 10 ** 3
v2 = 30 * 10 ** 3

# Wave functions
wave1 = phi_squared(x, sigma, t1, m_helium, v1)
wave2 = phi_squared(x, sigma, t1, m_helium, v2)

# Plot the wave functions
plt.plot(x, wave1, label='v1')
plt.plot(x, wave2, label='v2')
plt.legend()
plt.show()

# # ==============
# # Task c
# # ==============
#
# I = quad(phi_squared_integral, t1, t2, args=(sigma, t1, m_helium, v))
