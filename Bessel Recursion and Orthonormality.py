"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Orthonormality and recursion relation for Bessel functions
Author: Chitrak Roychowdhury
"""


import numpy as np
from scipy.special import jn,jvp
from scipy.misc import derivative
import scipy.integrate as sci

# Providing the value of m, start, stop, Np 
m = 20
start = 1e-2
stop = 10
Np = 1000 
z = np.linspace(start, stop, Np);  

# Logical case switch for different recursions to choose from
recb1=1; recb2=1; recb3=1; recb4=1; recb5=1; recb6=1; recb7=1;
print ('Compare maximum of |lhs-rhs| (L1 norm) to zero for m = ', m)

if(recb1): # z*J(m)'(z) = z*J(m-1)-m*J(m)
    lhs = z*jvp(m,z,1)
    rhs = z*jn(m-1,z)-m*jn(m,z)
    print ('Maximum of z*dJ(m)(z)/dz-z*J(m-1)+m*J(m) = ', abs(max(lhs-rhs)))

if(recb2): # 2*J(m)'(z) = J(m-1)(z)-J(m+1)(z)
    lhs = 2*jvp(m,z,1)
    rhs = jn(m-1,z)-jn(m+1,z)
    print ('Maximum of 2*dJ(m)(z)/dz-J(m-1)(z)+J(m+1)(z) = ', abs(max(lhs-rhs)))

if(recb3): # (2*m/z)*J(m)(z) = J(m+1)(z)+J(m-1)(z)
    lhs = np.divide(2*m*jn(m,z),z)
    rhs = jn(m+1,z)+jn(m-1,z)
    print ('Maximum of (2*m/z)*J(m)(z)-J(m+1)(z)+J(m-1)(z) = ', abs(max(lhs-rhs)))

if(recb4): # d(z^n*Jn(z))/dz = z^n*J(n-1)(z)
    def f(z):
        return jn(m,z)*(z**m)
    lhs = derivative(f, z, 1e-16)
    rhs = (z**m)*jn(m-1,z)
    print ('Maximum of d(z^n*J(n)(z))/dz-z^n*J(n-1)(z) = ', abs(max(lhs-rhs)))

if(recb5): # (z^(-m)*J(m)(z))' = -z^(-m)*J(m+1)(z)
    def f(z):
        return jn(m,z)*pow(z,-m)
    lhs = derivative(f, z, 1e-6)
    rhs = -pow(z,-m)*jn(m+1,z)
    print ('Maximum of d(z^(-m)*J(m)(z))/dz+z^(-m)*J(m+1)(z) = ', abs(max(lhs-rhs)))
    
if(recb6): # z*J(m)'(z) = m*J(m)-z*J(m+1)
    lhs = z*jvp(m,z,1)
    rhs = m*jn(m,z)-z*jn(m+1,z)
    print('Maximum of z*dJ(m)(z)/dz-m*J(m)+z*J(m+1) = ', abs(max(lhs-rhs)))


start = 0
stop = 1
Np = 1000 
z = np.linspace(start, stop, Np); 
a=int(input("Enter the value of a: "))
b=int(input("Enter the value of b: "))

if(recb7): # \int z*Jn(az)*Jn(bz) = { 0 if a!=b & 1/2*[J(n+1)(a)]^2 if a=b}
    I = sci.simps(z*jn(a*z,z)*jn(b*z,z),z)
    print ('Orthonormality :', I)
       

