"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Orthonormality and recursion relation for Chebyshev functions
Author: Chitrak Roychowdhury
"""


import numpy as np
from scipy.special import chebyt
from scipy.misc import derivative
import scipy.integrate as sci

# input the value of n, m, lower and upper limit
n=int(input("Enter the value of n:"))
m=int(input("Enter the value of m:"))
start = -1
stop = 1
Np= 1000
x = np.linspace(start, stop, Np)

# Create Poly1D Legendre polynomial and derivatives
tn = chebyt(n)
tm = chebyt(m)
tn1 = chebyt(n-1)
tn2 = chebyt(n+1)
tnprime = derivative(tn, x, 1e-6)  # spacing=10^-6
tn1prime = derivative(tn1, x, 1e-6)
tn2prime = derivative(tn2, x, 1e-6)

# Logical case switch for different recursion relations to choice from
rect1=1; rect2=1; rect3=1; rect4=1; 
print ('Compare maximum of |lhs-rhs| (L1 norm) to zero for n = ', n)

if(rect1): #T(n+1)(x) = 2*x*T(n)(x) - T(n-1)(x)
    lhs = tn2
    rhs = 2*x*tn(x)-(n-1)*tn1(x)
    print ('Maximum of T(n+1)(x)=2xT(n)(x)-T(n-1)(x) = ', abs(max(lhs-rhs)))

if(rect2): #(1-x^2)*T(n)'(x) = -n*x*T(n)(x) + n*T(n-1)(x)
    lhs = (1-x**2)*tnprime
    rhs = -n*x*tn(x)-n*tn1(x)
    print ('Maximum of (1-x^2)dT(n)(x)/dx=-nxT(n)(x)+nT(n-1)(x)  = ', abs(max(lhs-rhs)))

if(rect3): #2*T(n)(x) = T(n+1)'(x) - 2*x*T(n)'(x) + T(n-1)'(x)
    lhs = 2*tn(x)
    rhs = tn2prime + 2*x*tnprime + tn1prime
    print ('Maximum of 2*T(n)(x)=dT(n+1)(x)/dx-2xdT(n)(x)/dx + dT(n-1)(x)/dx  = ', abs(max(lhs-rhs)))

if(rect4): #\int T(n)(x)*T(m)(x)*(1-x**2)^(-1/2) = {0 if m!=n  pi/2 if m=n!=0  pi if m=n=0}
    I = sci.simps(tn(x)*tm(x),x)*(1-x**2)**(-0.5)
    print ('Orthonormality : = ', I)
    
    
        
