#!python
# -*- coding: utf-8 -*-

"""
...information about the code...
"""

__author__ = 'Astrid Moum'
__email__ = 'astridmo@nmbu.no'

import numpy as np
import math

a = 5 * 10**(-5)  # m
V0 = 5  # eV
m_elec = 0.109 * 10**(-31)  # kg
h_bar = 6.582 * 10**(-16)  # eV*s

E = np.arange(0.1, 1, 0.001)

k = np.sqrt(2*m_elec*E)/h_bar
kappa = np.sqrt(2*m_elec*(V0-E))/h_bar

def T(k, kappa, a):
    """Transmission"""
    return (1 + (((k**2 + kappa**2)**2) * (np.sinh(kappa*a)**2))/(4*(k**2)*(kappa**2)))


print(T(k, kappa, a))