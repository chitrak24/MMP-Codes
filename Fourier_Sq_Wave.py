"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Square Wave
Author: Chitrak Roychowdhury
"""

# Fourier series analysis for a sqaure wave function
# Using scipy.signal

import numpy as np
from scipy.signal import square
import matplotlib.pyplot as plt
from scipy.integrate import simps

period = 10 # Periodicity of the periodic function f(x)
freq = 5   # No of waves in time period L
Np = 100000
terms=100

# Generation of square wave
x=np.linspace(0,period,Np)
y=square(2.0*np.pi*x*freq/period,duty=0.5)

# Calculation of Fourier coefficients
a0=2./period*simps(y,x)
a_n=lambda m:2.0/period*simps(y*np.cos(2.*np.pi*m*x/period),x)
b_n=lambda m:2.0/period*simps(y*np.sin(2.*np.pi*m*x/period),x)

# sum of the series
s=a0/2.+sum([a_n(n)*np.cos(2.*np.pi*n*x/period)+b_n(n)*np.sin(2.*np.pi*n*x/period) for n in range(1,terms+1)])

#plotting
plt.plot(x,s, c='b', label="Fourier analysis for square wave")
plt.plot(x,y, c='r', label="Original square wave")
plt.xlabel("$x$")
plt.ylabel("$y=f(x)$")
plt.xlabel("x $\longrightarrow$")
plt.ylabel(r'$f(x)=\sum_{n=odd}^\infty\;\frac{4}{n\pi} sin(n\omega x)$')
plt.xticks()
plt.yticks()
plt.grid(True)
plt.xlim([0, period])
plt.axhline( color="green", linestyle="--")
plt.legend(loc='best',prop={'size':10})
plt.title("Sqaure wave signal analysis by Fouries series")
plt.savefig("fs_square.pdf")
plt.show()


# Fourier series analysis for a sqaure wave function
#Using User defined Function

A         = 4   # amplitude 
period    = 30  # periodicity
harmonics = 5   # Number of Harmonics
t = np.linspace(0,4*period,1000) # x-grid

# generate square waveform
def sqwave(t, period):
    return A*np.sign(np.sin(2*np.pi*t/period))
    
# fourier coefficients; an=0; bn=4/n*pi for odd-n.
def bn(n):                       
    if (n%2 != 0):
        return 10/(np.pi*n)
    else:
        return 0
    
# generate angular frequency
def wn(n, period):
    return (2*np.pi*n)/period

# Fourier series
def fourierSq(harmonics,t,period):
    summation = 0
    for i in range(1, harmonics):
        summation += A*bn(i)*np.sin(wn(i,period)*t)
    return summation

# Main 
y = []; f1 = []; f2 = []; f3 = []
for i in t:
    y.append(sqwave(i, period))
    f1.append(fourierSq(  harmonics,i,period))
    f2.append(fourierSq(4*harmonics,i,period))
    f3.append(fourierSq(8*harmonics,i,period))

# Plot
sg = A*square(2*np.pi*t/period)
plt.plot(t, sg, '-',  lw='4', color="red",    label="Original Square Waveform")
plt.plot(t,  y, '-o', lw='4', color="teal",   label=r'Signal$(A sgn(sin(\frac{2\pi t}{period}))$')
plt.plot(t, f1, '-*', lw='.5',color="blue",label=str(harmonics)+" harmonics")
plt.plot(t, f2, '-+', lw='.5',color="yellow",   label=str(4*harmonics)+" harmonics")
plt.plot(t, f3, '-x', lw='.5',color="olive",  label=str(8*harmonics)+" harmonics")
plt.title("Square Wave Fourier Series")
plt.legend(loc='best', prop={'size':16})
plt.xlabel('t $\longrightarrow$')
plt.xticks(size=14)
plt.ylabel(r'$f(t)=\sum_{n=odd}^\infty\;\frac{4A}{n\pi} sin(n\omega t)$', size=20)
plt.yticks(size=14)
plt.grid()
plt.savefig('fouriersq.pdf')
plt.show()
