"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Orthonormality and recursion relation for Laguerre functions
Author: Chitrak Roychowdhury
"""


import numpy as np
from scipy.special import laguerre
from scipy.misc import derivative
from scipy.integrate import quad
from scipy.integrate import simps

# Providing the value of m, start, stop, Np
m = 3
n = 2
start = 0
stop = 100
Np = 1000
x = np.linspace(start, stop, Np)

# Creating poly-1D Laguerre polynomials and their derivatives
ln = laguerre(n)
lm = laguerre(m)
lnp = laguerre(n+1)
lnq = laguerre(n-1)
lnprime = derivative(ln, x, 1e-6)
lnpprime = derivative(lnp, x, 1e-6)
lnqprime = derivative(lnq, x, 1e-6)

# Logical case switch for different recursions to choose from
recl1=1; recl2=1; recl3=1; recl4=1; recl5=1; recl6=1;
print ('Compare maximum of |lhs-rhs| (L1 norm) to zero for n = ', n)

if(recl1): # L'(n+1)=L'(n)-L(n)
    lhs = lnpprime
    rhs = lnprime - ln(x)
    print ('Maximum of dL(n+1)/dx-dL(n)/dx+L(n) = ', abs(max(lhs-rhs)))

if(recl2): # L(n)=e^x*L(n-1)(-x)
    lhs = ln(x)
    rhs = np.exp(x)*lnq(-x)
    print ('Maximum of L(n)=e^x*L(n-1)(-x) = ', abs(max(lhs-rhs)))

if(recl3): # (n+1)*L(n+1)=(2*n+1-x)*L(n)-n*L(n-1)
    lhs = (n+1)*lnp(x)
    rhs = (2*n+1-x)*ln(x)-n*lnq(x)
    print ('Maximum of (n+1)*L(n+1)=(2*n+1-x)*L(n)-n*L(n-1) = ', abs(max(lhs-rhs)))

if(recl4): # x*L'(n)=n*L(n)-n**2*L(n-1)
    lhs = x*lnprime
    rhs = n*ln(x) - n*lnq(x)
    print ('Maximum of x*dL(n)/dx=n*L(n)-n^2*L(n-1)  = ', abs(max(lhs-rhs)))

if(recl5): # L(n+1)=2*L(n)-L(n-1)-[1/(n+1)]*[(1+x)*L(n)-L(n-1)]
    lhs = lnp(x)
    rhs = 2*ln(x)-lnq(x)-(1/(n+1))*((1+x)*ln(x)-lnq(x))
    print ('Maximum of L(n+1)=2*L(n)-L(n-1)-[1/(n+1)]*[(1+x)*L(n)-L(n-1)]  = ', abs(max(lhs-rhs)))

#Checking Orthonormality:  # \int e^(-x)*L(n)(x)*L(m)(x) = (n!**2)*delta(nm)

from scipy.special import laguerre as l    
N=int(input("Enter degree of laguerre polynomial for checking orthonormality :"))

if(recl6):
    print ("Orthonormality: ")
    print ('m\t n\t value\t\t\t error\t ')
    print ('----------------------------------------------------------------------------')

    def f(x,M,N):
        return np.exp(-x)*l(M)(x)*l(N)(x)
    
    
    for M in range(0,N+1):
        I,err = quad(f,0,np.inf,args=(M,N))
        print ('%d\t%d\t%f\t\t%f '%(M,N,I,err))
        
       
        
----------------------------------------------------------------        
OUTPUT:        
----------------------------------------------------------------        
         
        
Compare maximum of |lhs-rhs| (L1 norm) to zero for n =  2
Maximum of dL(n+1)/dx-dL(n)/dx+L(n) =  4.774670560436789e-05
Maximum of L(n)=e^x*L(n-1)(-x) =  0.0
Maximum of (n+1)*L(n+1)=(2*n+1-x)*L(n)-n*L(n-1) =  1.1641532182693481e-10
Maximum of x*dL(n)/dx=n*L(n)-n^2*L(n-1)  =  5.309300104272552e-05
Maximum of L(n+1)=2*L(n)-L(n-1)-[1/(n+1)]*[(1+x)*L(n)-L(n-1)]  =  2.9103830456733704e-11
Enter degree of laguerre polynomial for checking orthonormality :6
Orthonormality: 
m	 n	 value			 error	 
----------------------------------------------------------------------------
0	6	0.000000		0.000000 
1	6	0.000000		0.000000 
2	6	0.000000		0.000000 
3	6	-0.000000		0.000000 
4	6	-0.000000		0.000000 
5	6	-0.000000		0.000000 
6	6	1.000000		0.000000 
