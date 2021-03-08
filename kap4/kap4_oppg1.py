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

# Setter opp uttrykk for høyre side
hs = np.sqrt(ksi0**2 - ksi**2)/ksi

# Venstre side
vs = np.tan(ksi)
print(vs)

# =========================
# Plotting
# # =======================
plt.figure(1)
plt.ylim(0, 25)
plt.xlim(0, max_value)
plt.plot(ksi, vs)
plt.plot(ksi, hs)
plt.show()




