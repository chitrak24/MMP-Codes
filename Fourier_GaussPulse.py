"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008
Description  : Fourier Series of Gaussian Pulse
Author       : Chitrak Roychowdhury
"""


import numpy as np
from scipy.signal import gausspulse
from scipy.integrate import simps
import matplotlib.pyplot as plt

L         = 1   # periodicity
harmonics = 30; # total number of term in series
x = np.linspace(-L,L,100) # x-grid


# Gauss Pulse
# Compute Fourier coefficients
a0 = 1.0/L*simps(gausspulse(x),x)
an = lambda n: 1.0/L*simps(gausspulse(x)*np.cos(n*np.pi*x/L),x)
bn = lambda n: 1.0/L*simps(gausspulse(x)*np.sin(n*np.pi*x/L),x)   

# Compute Fourier series
summgp = a0/2 + sum([an(k)*np.cos(k*np.pi*x/L) + bn(k)*np.sin(k*np.pi*x/L) for k in range(1,harmonics)])

#Plot
plt.plot(x, gausspulse(x),  color='cyan', ls='--', label = 'Signal')
plt.plot(x, summgp,   color='blue',       ls='dashdot', label = 'Computed')
plt.legend(loc='best')
plt.title("Gaussian Pulse", size=16);
plt.xlabel('$ x\longrightarrow $', size=16)
plt.xticks(size=14)
plt.ylabel('$ f(x)\longrightarrow $', size=20)
plt.yticks(size=14)
plt.grid(True)
plt.savefig("gausspulse.pdf")
plt.show()
