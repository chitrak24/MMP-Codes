"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Triangular Wave
Author: Chitrak Roychowdhury
"""

# Fourier series analysis for a Triangular wave function
# Using scipy.signal

import numpy as np
from scipy.signal import triang,sawtooth
import matplotlib.pyplot as plt
from scipy.integrate import simps

period = 10 # Periodicity of the periodic function f(x)
freq = 5   # No of waves in time period 
Np = 100000
terms=100

# Generation of Triangular wave
x=np.linspace(0,period,Np)
y=triang(Np)

# Calculation of Fourier coefficients
a0=2./period*simps(y,x)
a_n=lambda m:2.0/period*simps(y*np.cos(2.*np.pi*m*x/period),x)
b_n=lambda m:2.0/period*simps(y*np.sin(2.*np.pi*m*x/period),x)

# sum of the series
s=a0/2.+sum([a_n(n)*np.cos(2.*np.pi*n*x/period)+b_n(n)*np.sin(2.*np.pi*n*x/period) for n in range(1,terms+1)])

#plotting
plt.plot(x,s, color='gold', label="Fourier analysis for triangular wave")
plt.plot(x,y, color='cyan', label="Original triangular wave")
plt.xlabel("$x$")
plt.ylabel("$y=f(x)$")
plt.xlabel("x $\longrightarrow$")
plt.ylabel(r'$f(t)=\sum_{n=odd}^\infty\; \frac{8A}{n^2\pi^2}(-1)^{(n-1)/2} sin(n\omega t)$')
plt.xticks()
plt.yticks()
plt.grid(True)
plt.xlim([0, period])
plt.axhline( color="green", linestyle="--")
plt.legend(loc='best',prop={'size':10})
plt.title("triangular wave signal analysis by Fouries series")
plt.savefig("fs_triangular.pdf")
plt.show()


# Fourier series analysis for a Triangular wave function
#Using User defined Function

A         = 4   # amplitude 
period    = np.pi  # periodicity
harmonics = 3   # Number of Harmonics
x = np.linspace(-3*period,3*period,250) # x-grid

# generate Triangular waveform
def trwave(x, period):
    return A*2*np.arcsin(np.sin(np.pi*x/period))/np.pi
    
# fourier coefficients; an=0; bn=4/n*pi for odd-n.
def bn(n):
    if (n%2 != 0):
        return 8*pow(-1,(n-1)/2)/pow(np.pi*n,2)
    else:
        return 0
    
# generate angular frequency
def wn(n, period):
    return (2*np.pi*n)/period

# Fourier series
def fourierTr(harmonics,x):
    summation = 0
    for i in range(1, harmonics):
        summation += A*bn(i)*np.sin(i*np.pi*x/period)
    return summation

# Main 
y = []; f1 = []; f2 = []; f3 = []
for i in x:
    y.append(trwave(i, period))
    f1.append(fourierTr(  harmonics,i))
    f2.append(fourierTr(4*harmonics,i))
    f3.append(fourierTr(8*harmonics,i))

# Plot
tr = A*sawtooth(np.pi*(x+period/2)/period,width=0.5)
plt.plot(x, tr, '-',  lw='4', color="gold",   label="Original Triangular Waveform")
plt.plot(x,  y, '-o', lw='4', color="black",   label=r'Signal$(\frac{2A}{\pi}sin^{-1}(sin(\frac{\pi x}{period})$')
plt.plot(x, f1, '-*', lw='.5',color="red",    label=str(harmonics)+" harmonics")
plt.plot(x, f2, '-+', lw='.5',color="magenta", label=str(4*harmonics)+" harmonics")
plt.plot(x, f3, '-x', lw='.5',color="teal",    label=str(8*harmonics)+" harmonics")
plt.title("Triangular Wave Fourier Series")
plt.legend(loc='best')
plt.xlabel('x $\longrightarrow$')
plt.xticks(size=14)
plt.ylabel(r'$f(x)=\sum_{n=odd}^\infty\; \frac{8A}{n^2\pi^2}(-1)^{(n-1)/2} sin(n\omega x)$')
plt.yticks(size=14)
plt.grid()
plt.savefig('fouriertr.pdf')
plt.show()
