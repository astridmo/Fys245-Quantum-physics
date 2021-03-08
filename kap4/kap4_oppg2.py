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



def T_E_less_than_V0(m_elec, V0, E, a):
    """Transmission from task 4.26 for when E < V0"""
    k = np.sqrt(2 * m_elec * E) / h_bar
    kappa = np.sqrt(2 * m_elec * (V0 - E)) / h_bar
    return (1 + (((k**2 + kappa**2)**2) * (np.sinh(kappa*a)**2))/(4*(k**2)*(kappa**2)))


def T_E_larger_than_v0(m_elec, V0, E, a):
    """Transmission from task 4.28 for when E > V0"""
    k = np.sqrt(2 * m_elec * E) / h_bar
    k0 = np.sqrt(2 * m_elec * (E-V0)) / h_bar
    return (1 + (((k**2 - k0**2)**2) * (np.sin(k0*a)**2))/(4*(k**2)*(k0**2)))



E = np.arange(0.1, 3, 0.001)
print(T_E_less_than_V0(m_elec, V0, E, a))

E2 = np.arange(5.1, 8, 0.1)
print(T_E_larger_than_v0(m_elec, V0, E2, a))