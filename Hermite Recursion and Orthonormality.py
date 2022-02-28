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
from scipy.integrate import quad

# Providing the value of m, start, stop, Np 
n = 4
m = 5
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

#Checking Orthonormality:  # \int e^(-x^2)*H(n)(x)*H(m)(x) = (2^n)*(n!)*sqrt(pi)*delta(nm)

from scipy.special import hermite as h    
N=int(input("Enter degree of hermite polynomial for checking orthonormality :"))

if(rech6):
    print ("Orthonormality: ")
    print ('m\t n\t value\t\t error\t ')
    print ('----------------------------------------------------------------------------')

    def f(x,M,N):
        return np.exp(-x**2)*h(N)(x)*h(M)(x)
    
    for M in range(0,N+1):
        I,err = quad(f,-np.inf,np.inf,args=(M,N))
        print ('%d\t%d\t%f\t\t%f '%(M,N,I,err))

        
        
        
------------------------------------------------------------
OUTPUT:        
------------------------------------------------------------


        
Compare maximum of |lhs-rhs| (L1 norm) to zero for n =  4
Maximum of dH(n)(x)/dx - 2*n*H(n-1)(x)=  2.547281496845244e-09
Maximum of H(n)(-x)-(-1)^n*H(n)(x)=  0.0
Maximum of d^2H(n)(x)/dx^2 - 2*x*dH(n)(x)/dx - 2*n*H(n)(x)=  0.006259404165589899
Maximum of H(n+1)(x)=2*x*H(n)(x)-2*n*H(n-1)(x)=  1.4210854715202004e-14
Maximum of dH(n)(x)/dx -2*x*H(n)(x)+H(n+1)(x) =  2.547281496845244e-09
Enter degree of hermite polynomial for checking orthonormality :5
Orthonormality: 
m	 n	 value		 error	 
----------------------------------------------------------------------------
0	5	0.000000		0.000000 
1	5	0.000000		0.000000 
2	5	0.000000		0.000000 
3	5	0.000000		0.000000 
4	5	0.000000		0.000000 
5	5	6806.222787		0.000020 
