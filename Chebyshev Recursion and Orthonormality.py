"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Orthonormality and recursion relation for Chebyshev functions
Author: Chitrak Roychowdhury
"""


import numpy as np
from scipy.special import chebyt
from scipy.misc import derivative
from scipy.integrate import quad


# input the value of n, m, lower and upper limit
n=int(input("Enter the value of n:"))
m=int(input("Enter the value of m:"))
eps = 1e-3
start = -1+eps
stop = 1-eps
Np= 10000
x = np.linspace(start, stop, Np)

# Create Poly1D Legendre polynomial and derivatives
tn = chebyt(n)
tm = chebyt(m)
tn1 = chebyt(n-1)
tn2 = chebyt(n+1)
tnprime = derivative(tn, x, 1e-5)  # spacing=10^-5
tn1prime = derivative(tn1, x, 1e-5)
tn2prime = derivative(tn2, x, 1e-5)

# Logical case switch for different recursion relations to choice from
rect1=1; rect2=1; rect3=1; rect4=1; 
print ('Compare maximum of |lhs-rhs| (L1 norm) to zero for n = ', n)

if(rect1): #T(n+1)(x) = 2*x*T(n)(x) - T(n-1)(x)
    lhs = tn2(x)
    rhs = 2*x*tn(x) - tn1(x)
    print ('Maximum of T(n+1)(x)=2xT(n)(x)-T(n-1)(x) = ', abs(max(lhs-rhs)))

if(rect2): #(1-x^2)*T(n)'(x) = -n*x*T(n)(x) + n*T(n-1)(x)
    lhs = (1-x**2)*tnprime
    rhs = -n*x*tn(x) + n*tn1(x)
    print ('Maximum of (1-x^2)dT(n)(x)/dx=-nxT(n)(x)+nT(n-1)(x)  = ', abs(max(lhs-rhs)))

if(rect3): #2*T(n)(x) = T(n+1)'(x) - 2*x*T(n)'(x) + T(n-1)'(x)
    lhs = 2*tn(x)
    rhs = tn2prime - 2*x*tnprime + tn1prime
    print ('Maximum of 2*T(n)(x)=dT(n+1)(x)/dx-2xdT(n)(x)/dx + dT(n-1)(x)/dx  = ', abs(max(lhs-rhs)))

#Checking Orthonormality

from scipy.special import chebyt as t    
N=int(input("Enter degree of chebyshev polynomial for checking orthonormality :"))

if(rect4):
    print ("Orthonormality: ")
    print ('m\t n\t delta_mn\t\t error\t ')
    print ('----------------------------------------------------------------------------')

    def f(x,M,N):
        return ((t(M)(x)*t(N)(x))/np.sqrt(1-x**2))
    
    for M in range(0,N+1):
        I,err = quad(f,-1,1,args=(M,N))
        print ('%d\t%d\t%f\t\t%f '%(M,N,I,err))

        
        
---------------------------------------------        
OUTPUT:        
---------------------------------------------
        
        
        
Enter the value of n:5
Enter the value of m:6
Compare maximum of |lhs-rhs| (L1 norm) to zero for n =  5
Maximum of T(n+1)(x)=2xT(n)(x)-T(n-1)(x) =  6.661338147750939e-16
Maximum of (1-x^2)dT(n)(x)/dx=-nxT(n)(x)+nT(n-1)(x)  =  3.0758635638505893e-09
Maximum of 2*T(n)(x)=dT(n+1)(x)/dx-2xdT(n)(x)/dx + dT(n-1)(x)/dx  =  1.9901232217733877e-08
Enter degree of chebyshev polynomial for checking orthonormality :9
Orthonormality: 
m	 n	 delta_mn		 error	 
----------------------------------------------------------------------------
0	9	0.000000		0.000000 
1	9	0.000000		0.000000 
2	9	0.000000		0.000000 
3	9	0.000000		0.000000 
4	9	0.000000		0.000000 
5	9	0.000000		0.000000 
6	9	0.000000		0.000000 
7	9	0.000000		0.000000 
8	9	0.000000		0.000000 
9	9	1.570796		0.000000 
