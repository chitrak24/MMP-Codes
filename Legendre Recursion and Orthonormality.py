"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Orthonormality and recursion relation for Legendre functions
Author: Chitrak Roychowdhury
"""


import numpy as np
from scipy.special import legendre
from scipy.misc import derivative
import scipy.integrate as sci

# input the value of n, m, lower and upper limit
n=int(input("Enter the value of n:"))
m=int(input("Enter the value of m:"))
start = -1
stop = 1
Np= 10000000
x = np.linspace(start, stop, Np)

# Create Poly1D Legendre polynomial and derivatives
pn = legendre(n)
pm = legendre(m)
pnm1 = legendre(n-1)
pnm2 = legendre(n-2)
pnp1 = legendre(n+1)
pnprime = derivative(pn, x, 1e-6)  # spacing=10^-6
pnm1prime = derivative(pnm1, x, 1e-6)
pnp1prime = derivative(pnp1, x, 1e-6)
pnm1prime = derivative(pnm1, x, 1e-6)

# Logical case switch for different recursion relations to choice from
recl1=1; recl2=1; recl3=1; recl4=1; recl5=1; recl6=1; recl7=1; recl8=1; recl9=1;
print ('Compare maximum of |lhs-rhs| (L1 norm) to zero for n = ', n)

if(recl1): # \int Pn(x)Pm(x) = 2/(2n+1)*delta(nm)
    I = sci.simps(pn(x)*pm(x),x)*(2.0*n+1)/2.0
    print ('Orthonormality : \int P_',n,'(x)P_',m,'(x)dx = ', I)

if(recl2): #nPn(x) = (2n-1)xP(n-1)(x) - (n-1)P(n-2)(x)
    lhs = n*pn(x)
    rhs = (2*n-1)*x*pnm1(x)-(n-1)*pnm2(x)
    print ('Maximum of nPn(x)-(2n-1)xPn(x)+(n-1)P(n-2)(x) = ', abs(max(lhs-rhs)))

if(recl3): #(n+1)P(n+1)(x) = (2n+1)xPn(x) - nP(n-1)(x)
    lhs = (n+1)*pnp1(x)
    rhs = (2*n+1)*x*pn(x)-n*pnm1(x)
    print ('Maximum of (n+1)P(n+1)(x)-(2n+1)xPn(x)+nP(n-1)(x) = ', abs(max(lhs-rhs)))

if(recl4): #(1-x^2)Pn'(x) = n(P(n-1)(x) - xPn(x))
    lhs = (1-x**2)*pnprime
    rhs = n*(pnm1(x)-x*pn(x))
    print ('Maximum of (1-x^2)dPn(x)/dx-n[P(n-1)(x)-xPn(x)] = ', abs(max(lhs-rhs)))

if(recl5): #nPn(x) = xPn'(x) - P(n-1)'(x)
    lhs = n*pn(x)
    rhs = x*pnprime-pnm1prime
    print ('Maximum of nPn(x)-xdPn(x)/dx+dP(n-1)(x)/dx = ', abs(max(lhs-rhs)))

if(recl6): #Pn'(x) = xP(n-1)'(x) + nP(n-1)(x)
    lhs = pnprime
    rhs = x*pnm1prime + n*pnm1(x)
    print ('Maximum of dPn(x)/dx-xdP(n-1)(x)/dx-nP(n-1)(x) = ', abs(max(lhs-rhs)))

if(recl7): #(2n+1)Pn(x) = P(n+1)'(x) - P(n-1)'(x)
    lhs = (2*n+1)*pn(x)
    rhs = pnp1prime - pnm1prime
    print ('Maximum of (2n+1)Pn(x)+dP(n-1)(x)/dx-dP(n+1)(x)/dx = ', abs(max(lhs-rhs)))

if(recl8): #(n+1)Pn(x) = P(n+1)'(x) - xP(n)'(x)
    lhs = (n+1)*pn(x)
    rhs = pnp1prime - x*pnprime
    print ('Maximum of (n+1)Pn(x)+xdP(n)(x)/dx-dP(n+1)(x)/dx = ', abs(max(lhs-rhs)))

if(recl9): #(n+1)[xP(n)(x) - P(n+1)(x)]=(1-x^2)P(n)'(x)
    lhs = (n+1)*x*pn(x) - (n+1)*pnp1(x)
    rhs = (1-x**2)*pnprime
    print ('Maximum of (n+1)[xP(n)(x)-P(n+1)(x)]-(1-x^2)dP(n)(x)/dx = ', abs(max(lhs-rhs)))
    
