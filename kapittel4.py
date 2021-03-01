#!python
# -*- coding: utf-8 -*-

"""
...information about the code...
"""

__author__ = 'Astrid Moum'
__email__ = 'astridmo@nmbu.no'

import math
import matplotlib.pyplot as plt

a = 2
h_bar = 1.05457 * 10 ^ (-34)
m_elec = 9.10938356 * 10 ^ (-31)

E = 10  # energy of ??

tsi = math.sqrt(2 * m_elec * E) * a * (1 / (h_bar * 2))
tsi0 = 10
# tsi0 = (a/2) * math.sqrt(2*m_elec*V0)/h_bar

tan_tsi = math.tan(tsi)


