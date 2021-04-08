# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 09:14:11 2021

@author: Johan
"""

# Impoerterer
import numpy as np
import matplotlib.pyplot as plt

# Definerer tidssteg
dt = 2.25 * 10 ** (-19)  # s
t_start = 0  # [s]
t_stop = 10 ** (-14)  # [s]
N_sim = int((t_stop - t_start) / dt) + 1  # Antall tidssteg (simulasjoner) som gjøres

# Definerer variabler
sigma = 1 * 10 ** (-8)  # m
h_bar = 1.055 * 10 ** (-34)  # eV s
m = 9.11 * 10 ** (-31)  # kg
E = 0.2 * 1.602 * 10 ** (-19)  # eV
L = 200 * 10 ** (-9)  # m

x0 = 50 * 10 ** (-9)  # m
dx = 0.15 * 10 ** (-9)  # m
x = np.arange(0, L, dx)  # m

# t0 = 0 #s
##dt = 2.25 * 10**(-19) #s
# t1 = t0 + 5000*dt

# t = np.arange(t0, t1, dt)


# Finner løsningen i t = 0
temp_p = (1 / (np.pi ** (1 / 4) * np.sqrt(sigma))) * np.exp(
    -((x - x0) ** 2 / (2 * sigma ** 2)) + 1j * (np.sqrt(2 * m * E) / h_bar) * (x - x0))
x_len = int(len(temp_p))

# Definerer arrays til å lagre verdier i
Phi_array = np.zeros((N_sim, x_len), dtype=complex)
Phi2_array = np.zeros((N_sim, x_len), dtype=complex)
time_array = np.zeros(N_sim)


# Definerer trinnpotensialet
def V(x):
    """
    Potensialfunksjonen: x == 0 og x == L settes til 0, da verdiene her settes til 0 uansett
    """

    if x > L / 2 and x < L:
        return (0.16 * 1.602 * 10 ** (-19))
    else:
        return 0


# def Phi_neste(step):
#
#    Phi_forrige = temp_p[step]
#    delta = ((dt/complex(0, h_bar)) * (temp_p[step]*V(x[step]) - (h_bar**2/(2*m*dx**2))*(temp_p[step+1]-2*temp_p[step]+temp_p[step-1])))
#    Phi_neste[step] = temp_p[step] + delta


# Regner for hvert tidssteg
for k in range(0, N_sim):
    t_k = k * dt

    # Forskjell mellom tidssteg
    delta = np.zeros(len(temp_p), dtype=complex)
    Phi_neste = np.zeros(len(temp_p), dtype=complex)

    for step in range(int(len(temp_p))):
        step -= 1
        if step <= 0:
            delta[step] = 0
        else:
            delta[step] = ((dt / complex(0, h_bar)) * (temp_p[step] * V(x[step]) - (h_bar ** 2 / (2 * m * dx ** 2)) * (
                        temp_p[step + 1] - 2 * temp_p[step] + temp_p[step - 1])))

    # Regner verdi i neste tidssteg og setter endepunktene til null
    Phi_neste = temp_p + delta
    Phi_neste[0] = 0
    Phi_neste[-1] = 0

    # Opphøyer phi i annen
    Phi2 = Phi_neste * np.conj(Phi_neste)

    # Oppdaterer arrayene
    time_array[k] = t_k
    Phi_array[k] = Phi_neste
    Phi2_array[k] = Phi2

    # Forbereder for neste loop
    temp_p = Phi_neste

# plt.figure(1)
# plt.plot(x, Phi_array[0])
# plt.plot(x, Phi_array[-1])

plt.figure(2)
plt.plot(x, Phi2_array[0])
plt.plot(x, Phi2_array[-1])
