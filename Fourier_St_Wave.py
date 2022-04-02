"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Sawtooth Wave
Author: Chitrak Roychowdhury
"""

# Fourier series analysis for a sawtooth wave function
# Using scipy.signal

import numpy as np
from scipy.signal import sawtooth
import matplotlib.pyplot as plt
from scipy.integrate import simps

period = 10 # Periodicity of the periodic function f(x)
freq = 5   # No of waves in time period L
Np = 100000
terms=100

# Generation of square wave
x=np.linspace(0,period,Np)
y=sawtooth(2.0*np.pi*x*freq/period)

# Calculation of Fourier coefficients
a0=2./period*simps(y,x)
a_n=lambda m:2.0/period*simps(y*np.cos(2.*np.pi*m*x/period),x)
b_n=lambda m:2.0/period*simps(y*np.sin(2.*np.pi*m*x/period),x)

# sum of the series
s=a0/2.+sum([a_n(n)*np.cos(2.*np.pi*n*x/period)+b_n(n)*np.sin(2.*np.pi*n*x/period) for n in range(1,terms+1)])

#plotting
plt.plot(x,s, color='cyan', label="Fourier analysis for sawtooth wave")
plt.plot(x,y, color='black', label="Original sawtooth wave")
plt.xlabel("$x$")
plt.ylabel("$y=f(x)$")
plt.xlabel("x $\longrightarrow$")
plt.ylabel(r'$f(x)=\sum_{n=1}^\infty\;\frac{2}{n\pi}(-1)^{n+1} sin(n\omega x)$')
plt.xticks()
plt.yticks()
plt.grid(True)
plt.xlim([0, period])
plt.axhline( color="green", linestyle="--")
plt.legend(loc='best',prop={'size':10})
plt.title("Sawtooth wave signal analysis by Fouries series")
plt.savefig("fs_sawtooth.pdf")
plt.show()


# Fourier series analysis for a Sawtooth wave function
#Using User defined Function

A         = 4   # amplitude 
period    = np.pi  # periodicity
harmonics = 3   # Number of Harmonics
t = np.linspace(-3*period,3*period,250) # x-grid

# generate square waveform
def stwave(t, period):
    return A*2*(t/period - np.floor(.5+t/period))
    
# fourier coefficients; an=0; bn=4/n*pi for odd-n.
def bn(n):                       
    return pow(-1,n+1)*2/(np.pi*n)
    
# generate angular frequency
def wn(n, period):
    return (2*np.pi*n)/period

# Fourier series
def fourierSt(harmonics,t,period):
    summation = 0
    for i in range(1, harmonics):
        summation += A*bn(i)*np.sin(wn(i,period)*t)
    return summation

# Main 
y = []; f1 = []; f2 = []; f3 = []
for i in t:
    y.append(stwave(i, period))
    f1.append(fourierSt(  harmonics,i,period))
    f2.append(fourierSt(4*harmonics,i,period))
    f3.append(fourierSt(8*harmonics,i,period))

# Plot
sg = A*sawtooth(2*(t-period/2))
plt.plot(t, sg, '-',  lw='4', color="olive",   label="Original Sawtooth Waveform")
plt.plot(t,  y, '-o', lw='4', color="black",   label=r'Signal$(2A(\frac{t}{period} - floor({\frac{1}{2}+\frac{t}{period})}$')
plt.plot(t, f1, '-*', lw='.5',color="cyan",    label=str(harmonics)+" harmonics")
plt.plot(t, f2, '-+', lw='.5',color="magenta", label=str(4*harmonics)+" harmonics")
plt.plot(t, f3, '-x', lw='.5',color="teal",    label=str(8*harmonics)+" harmonics")
plt.title("Sawtooth Wave Fourier Series")
plt.legend(loc='best', prop={'size':16})
plt.xlabel('t $\longrightarrow$')
plt.xticks(size=14)
plt.ylabel(r'$f(t)=\sum_{n=1}^\infty\;\frac{2}{n\pi}(-1)^{n+1} sin(n\omega t)$', size=20)
plt.yticks(size=14)
plt.grid()
plt.savefig('fourierst.pdf')
plt.show()
