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
from scipy.integrate import simps

# input the value of n, m, lower and upper limit
n=int(input("Enter the value of n:"))
m=int(input("Enter the value of m:"))
start = -1
stop = 1
Np= 100
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

if(recl1): #nPn(x) = (2n-1)xP(n-1)(x) - (n-1)P(n-2)(x)
    lhs = n*pn(x)
    rhs = (2*n-1)*x*pnm1(x)-(n-1)*pnm2(x)
    print ('Maximum of nPn(x)-(2n-1)xPn(x)+(n-1)P(n-2)(x) = ', abs(max(lhs-rhs)))

if(recl2): #(n+1)P(n+1)(x) = (2n+1)xPn(x) - nP(n-1)(x)
    lhs = (n+1)*pnp1(x)
    rhs = (2*n+1)*x*pn(x)-n*pnm1(x)
    print ('Maximum of (n+1)P(n+1)(x)-(2n+1)xPn(x)+nP(n-1)(x) = ', abs(max(lhs-rhs)))

if(recl3): #(1-x^2)Pn'(x) = n(P(n-1)(x) - xPn(x))
    lhs = (1-x**2)*pnprime
    rhs = n*(pnm1(x)-x*pn(x))
    print ('Maximum of (1-x^2)dPn(x)/dx-n[P(n-1)(x)-xPn(x)] = ', abs(max(lhs-rhs)))

if(recl4): #nPn(x) = xPn'(x) - P(n-1)'(x)
    lhs = n*pn(x)
    rhs = x*pnprime-pnm1prime
    print ('Maximum of nPn(x)-xdPn(x)/dx+dP(n-1)(x)/dx = ', abs(max(lhs-rhs)))

if(recl5): #Pn'(x) = xP(n-1)'(x) + nP(n-1)(x)
    lhs = pnprime
    rhs = x*pnm1prime + n*pnm1(x)
    print ('Maximum of dPn(x)/dx-xdP(n-1)(x)/dx-nP(n-1)(x) = ', abs(max(lhs-rhs)))

if(recl6): #(2n+1)Pn(x) = P(n+1)'(x) - P(n-1)'(x)
    lhs = (2*n+1)*pn(x)
    rhs = pnp1prime - pnm1prime
    print ('Maximum of (2n+1)Pn(x)+dP(n-1)(x)/dx-dP(n+1)(x)/dx = ', abs(max(lhs-rhs)))

if(recl7): #(n+1)Pn(x) = P(n+1)'(x) - xP(n)'(x)
    lhs = (n+1)*pn(x)
    rhs = pnp1prime - x*pnprime
    print ('Maximum of (n+1)Pn(x)+xdP(n)(x)/dx-dP(n+1)(x)/dx = ', abs(max(lhs-rhs)))

if(recl8): #(n+1)[xP(n)(x) - P(n+1)(x)]=(1-x^2)P(n)'(x)
    lhs = (n+1)*x*pn(x) - (n+1)*pnp1(x)
    rhs = (1-x**2)*pnprime
    print ('Maximum of (n+1)[xP(n)(x)-P(n+1)(x)]-(1-x^2)dP(n)(x)/dx = ', abs(max(lhs-rhs)))
    

#Checking Orthonormality:  # \int e^(-x^2)*H(n)(x)*H(m)(x) = (2^n)*(n!)*sqrt(pi)*delta(nm)

from scipy.special import legendre as p    
N=int(input("Enter degree of legendre polynomial for checking orthonormality :"))

if(recl9):
    print ("Orthonormality: ")
    print ('m\t n\t\t delta_mn\t ')
    print ('----------------------------------------------------------------------------')
    
    for M in range(0,N+2):
        f = p(M)(x)*p(N)(x)
        result = ((2*N+1)/2)*simps(f,x)
        print ('%d\t%d\t\t%d '%(M,N,result))

        
     
    
------------------------------------------------------------    
OUTPUT:
------------------------------------------------------------    
    
    
Enter the value of n:5
Enter the value of m:6
Compare maximum of |lhs-rhs| (L1 norm) to zero for n =  5
Maximum of nPn(x)-(2n-1)xPn(x)+(n-1)P(n-2)(x) =  1.7763568394002505e-15
Maximum of (n+1)P(n+1)(x)-(2n+1)xPn(x)+nP(n-1)(x) =  2.220446049250313e-15
Maximum of (1-x^2)dPn(x)/dx-n[P(n-1)(x)-xPn(x)] =  1.4847501006443053e-10
Maximum of nPn(x)-xdPn(x)/dx+dP(n-1)(x)/dx =  6.988898348936345e-10
Maximum of dPn(x)/dx-xdP(n-1)(x)/dx-nP(n-1)(x) =  6.988898348936345e-10
Maximum of (2n+1)Pn(x)+dP(n-1)(x)/dx-dP(n+1)(x)/dx =  1.537557636765996e-09
Maximum of (n+1)Pn(x)+xdP(n)(x)/dx-dP(n+1)(x)/dx =  8.386678018723615e-10
Maximum of (n+1)[xP(n)(x)-P(n+1)(x)]-(1-x^2)dP(n)(x)/dx =  2.3881585597962385e-10
Enter degree of legendre polynomial for checking orthonormality :10
Orthonormality: 
m	 n		 delta_mn	 
----------------------------------------------------------------------------
0	10		0 
1	10		0 
2	10		0 
3	10		0 
4	10		0 
5	10		0 
6	10		0 
7	10		0 
8	10		0 
9	10		0 
10	10		1 
11	10		0 
   
    
    
    
    
    
 
