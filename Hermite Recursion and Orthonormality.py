"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Orthonormality and recursion relation for Hermite functions
Author: Chitrak Roychowdhury
"""


import numpy as np
import math
from scipy.special import hermite
from scipy.misc import derivative
import scipy.integrate as sci

# Providing the value of m, start, stop, Np 
n = 4
m = 3
start = -0.5
stop = 0.5
Np = 1000
x = np.linspace(start, stop, Np)

#creating poly-1D Hermite polynomials and their derivatives
hn = hermite(n)
hm = hermite(m)
hnp = hermite(n-1)
hnq = hermite(n+1)
hnprime = derivative(hn, x, 1e-6)
hndprime = derivative(hn, x, 1e-6, 2)

# Logical case switch for different recursions to choose from
rech1=1; rech2=1; rech3=1; rech4=1; rech5=1; rech6=1;
print ('Compare maximum of |lhs-rhs| (L1 norm) to zero for n = ', n)

if(rech1): # H(n)'(x)=2*n*H(n-1)(x)
    lhs = hnprime
    rhs = 2*n*hnp(x)
    print('Maximum of dH(n)(x)/dx - 2*n*H(n-1)(x)= ', abs(max(lhs-rhs)))

if(rech2): # H(n)(-x)=(-1)^n*H(n)(x)
    lhs = hn(-x)
    rhs = pow(-1,n)*hn(x)
    print('Maximum of H(n)(-x)-(-1)^n*H(n)(x)= ', abs(max(lhs-rhs)))

if(rech3): # H(n)''(x)=2*x*H(n)'(x)-2*n*H(n)(x)
    lhs = hndprime
    rhs = 2*x*hnprime-2*n*hn(x)
    print('Maximum of d^2H(n)(x)/dx^2 - 2*x*dH(n)(x)/dx - 2*n*H(n)(x)= ', abs(max(lhs-rhs)))

if(rech4): # H(n+1)(x)=2*x*H(n)(x)-2*n*H(n-1)(x)
    lhs = hnq(x)
    rhs = 2*x*hn(x)-2*n*hnp(x)
    print('Maximum of H(n+1)(x)=2*x*H(n)(x)-2*n*H(n-1)(x)= ', abs(max(lhs-rhs)))

if(rech5): # H(n)'(x)=2*x*H(n)(x)-H(n+1)(x)
    lhs = hnprime
    rhs = 2*x*hn(x)-hnq(x)
    print('Maximum of dH(n)(x)/dx -2*x*H(n)(x)+H(n+1)(x) = ', abs(max(lhs-rhs)))

if(rech6): # \int e^(-x^2)*H(n)(x)*H(m)(x) = (2^n)*(n!)*sqrt(pi)*delta(nm)
    I = sci.simps(np.exp(-x**2)*hn(x)*hm(x),x)*np.divide(1,pow(2,n)*math.factorial(n)*np.sqrt(np.pi))
    print ('Orthonormality : = ', I)
  
