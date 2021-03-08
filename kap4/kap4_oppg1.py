#!python
# -*- coding: utf-8 -*-

"""
...information about the code...
"""

__author__ = 'Astrid Moum'
__email__ = 'astridmo@nmbu.no'

import math
import matplotlib.pyplot as plt
import numpy as np

# Definerer en verdi for ksi0
v = 5 # v = 1 gir 1 krysningspunkt, v = 2 gir 2 krysningspunkter osv.


# Definerer en array som definerer ksi
max_value = (v-0.5)*math.pi   # Max value of ksi
ksi = np.linspace(0, max_value, 1000)  # Definerer ksi
ksi0 = math.pi

#print(ksi)
# Setter opp uttrykk for høyre side
try:
    hs = np.sqrt(ksi0**2 - ksi**2)/ksi
except RuntimeWarning:
    pass

# Venstre side
try:
    vs = np.tan(ksi)
except RuntimeWarning:
    pass
print(vs)

vs = vs[:-1]

#plt.plot(ksi[:-1], vs)
plt.plot(ksi, hs)
plt.show()




