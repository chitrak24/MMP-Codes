"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Orthonormality and recursion relation for Laguerre functions
Author: Chitrak Roychowdhury
"""


import numpy as np
from scipy.special import laguerre
from scipy.misc import derivative
import scipy.integrate as sci

# Providing the value of m, start, stop, Np
m = 3
n = 2
start = 0
stop = 5.0
Np = 1000
x = np.linspace(start, stop, Np)

# Creating poly-1D Laguerre polynomials and their derivatives
ln = laguerre(n)
lm = laguerre(m)
lnp = laguerre(n+1)
lnq = laguerre(n-1)
lnprime = derivative(ln, x, 1e-6)
lnqprime = derivative(lnq, x, 1e-6)

# Logical case switch for different recursions to choose from
recl1=1; recl2=1; recl3=1; recl4=1; recl5=1; recl6=1;
print ('Compare maximum of |lhs-rhs| (L1 norm) to zero for n = ', n)

if(recl1): # L'(n)=n*L'(n-1)-n*L(n-1)
    lhs = lnprime
    rhs = n*lnqprime -n*lnq(x)
    print ('Maximum of dL(n)/dx-n*dL(n-1)/dx+n*L(n-1) = ', abs(max(lhs-rhs)))

if(recl2): # L(n+1)+(n^2)*L(n-1)=-(x-2*n-1)*L(n)
    lhs = lnp(x) + (n**2)*lnq(x)
    rhs = -(x-2*n-1)*ln(x)
    print ('Maximum of L(n+1)+(n^2)*L(n-1)+(x-2*n-1)*L(n) = ', abs(max(lhs-rhs)))

if(recl3): # (n+1)*L(n+1)=(2*n+1-x)*L(n)-n*L(n-1)
    lhs = (n+1)*lnp(x)
    rhs = (2*n+1-x)*ln(x)-n*lnq(x)
    print ('Maximum of (n+1)*L(n+1)=(2*n+1-x)*L(n)-n*L(n-1) = ', abs(max(lhs-rhs)))

if(recl4): # x*L'(n)=n*L(n)-n*L(n-1)
    lhs = x*lnprime
    rhs = n*ln(x) - n*lnq(x)
    print ('Maximum of x*dL(n)/dx=n*L(n)-n*L(n-1)  = ', abs(max(lhs-rhs)))

if(recl5): # L(n+1)=2*L(n)-L(n-1)-[1/(n-1)]*[(1+x)*L(n)-L(n-1)]
    lhs = lnp(x)
    rhs = 2*ln(x)-lnq(x)-(1/(n-1))*((1+x)*ln(x)-lnq(x))
    print ('Maximum of L(n+1)=2*L(n)-L(n-1)-[1/(n-1)]*[(1+x)*L(n)-L(n-1)]  = ', abs(max(lhs-rhs)))

if(recl6): # \int e^(-x)*L(n)(x)*L(m)(x) = (n!**2)*delta(nm)
    I = sci.simps(np.exp(-x)*lm(x)*ln(x),x)
    print ('Orthonormality : ', I)
