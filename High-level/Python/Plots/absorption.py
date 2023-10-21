"""
Atlantic Ocean seawater,  (T = 10°C, S = 35 ppt, pH = 7.8, D=3)
"""
import matplotlib.pyplot as plt
import numpy as np

TEMP = 10
ABS_TEMP = TEMP + 273 # Used as Theta symbol in the equations. 
PRESSURE = 101325
SALINITY = 35 # ppt: parts per thousands
DEPTH = 3
pH = 7.8

sound_speed = 1412 + 3.21*TEMP + 1.19*SALINITY + 0.167 * DEPTH

freq_kHz = np.logspace(-2, 3, num=1000)

def viscous_absorption(freq_kHz): 
    if TEMP <= 20:
        A_3 = 4.937 * pow(10, -4) - 2.59 * pow(10, -5) * TEMP + 9.11 * pow (10, -7) * TEMP ** 2 - 1.50* pow(10, -8)*TEMP ** 3
        # What about P_3 when temp <= 20? Could not find that in the book
        P_3 = 1
    else: 
        A_3 = 3.964 * pow(10, -4) - 1.146 * pow(10, -5) * TEMP + 1.45 * pow (10, -7) * TEMP ** 2 - 6.5* pow(10, -10)*TEMP ** 3
        P_3 = 1 - 3.83*DEPTH * pow(10, -5) + DEPTH ** 2 * pow(10, -10)

    return A_3 * P_3 * pow(freq_kHz, 2)

# Boric acid contribution
A_1 = (8.86 / sound_speed) * pow(10, 0.78*pH-5)
P_1 = 1
f_1 = 2.8 * pow(SALINITY/35, 0.5) * pow(10, 4-1245/ABS_TEMP)

# Magnesium sulfate contribution
A_2 = 21.44*(SALINITY/sound_speed)*(1 + 0.025*TEMP)
P_2 = 1 - 1.37*DEPTH * pow(10, -4) + 6.2 * DEPTH ** 2 * pow(10, -9)
f_2 = (8.17 * pow(10, 8-1990/ABS_TEMP))/ (1 + 0.0018 * (SALINITY - 35))


# Boric acid contribution
BoricAcid = (A_1 * P_1 * f_1 * freq_kHz ** 2) / (freq_kHz ** 2 + f_1 ** 2)

# Magnesium sulfate contribution
MagSulf = (A_2 * P_2 * f_2 * freq_kHz ** 2) / (freq_kHz ** 2 + f_2 ** 2)

# Viscous absoroption 3rd terms in equation
freshViscous = viscous_absorption(freq_kHz)

# The combined equation for all the 3 terms in equation (3.27)
alphaTot = BoricAcid + MagSulf +  freshViscous

plt.figure(figsize=(10, 6))

plt.plot(freq_kHz, alphaTot, label='Salt Water', color='black', linestyle='solid')
plt.plot(freq_kHz, freshViscous, label='Fresh water ', color='blue', linestyle='--', linewidth=1)
plt.plot(freq_kHz, MagSulf, label=r'MgSO$_4$', color='green', linestyle=':', linewidth=1)
plt.plot(freq_kHz, BoricAcid, label=r'B(OH)$_3$', color='red', linestyle='-.', linewidth=1)


plt.xlabel('Frequency (kHz)')
plt.ylabel('Attenuation - dB/km')

plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-2, 1e3)
plt.ylim(1e-6, 1e4)

plt.grid(which='both', linestyle='--', linewidth=0.5)

plt.legend()

# plt.savefig('absorption.pdf', format='pdf')
plt.show()