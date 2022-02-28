"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Orthonormality and recursion relation for Bessel functions
Author: Chitrak Roychowdhury
"""


import numpy as np
from scipy.special import jn,jvp,jn_zeros
from scipy.misc import derivative
import scipy.integrate as sci
import random

# Providing the value of m, start, stop, Np 
m = 20
start = 1e-2
stop = 10
Np = 1000
z = np.linspace(start, stop, Np);  

# Logical case switch for different recursions to choose from
recb1=1; recb2=1; recb3=1; recb4=1; recb5=1; recb6=1;
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


#Checking Orthonormality:  # \int z*Jn(az)*Jn(bz) = { 0 if a!=b & 1/2*[J(n+1)(a)]^2 if a=b}

n=int(input('Enter the degree of the polynomial:'))
nor=int(input('Enter the total number of roots:'))
a=jn_zeros(n,nor)                                    # Calculation of roots of the Jn(x)=0
z=np.linspace(0,1,100)                               # limit of integration

print ('---------------------------------------------------------')
print (' a=\t\t\t',' b=\t\t',' \tI=')
print ('---------------------------------------------------------')

i=random.choice(a)                                   #selecting a random root

for j in a:
 fun=z*jn(n,i*z)*jn(n,j*z)
 rhs=0.5*(jn(n+1,i))**2
 I=sci.simps(fun,z)/rhs
 print('%f\t\t %f\t\t %f'%(i,j,I))
 



-------------------------------------------------------
OUTPUT:
-------------------------------------------------------

Compare maximum of |lhs-rhs| (L1 norm) to zero for m =  20
Maximum of z*dJ(m)(z)/dz-z*J(m-1)+m*J(m) =  1.0164395367051604e-19
Maximum of 2*dJ(m)(z)/dz-J(m-1)(z)+J(m+1)(z) =  0.0
Maximum of (2*m/z)*J(m)(z)-J(m+1)(z)+J(m-1)(z) =  2.371692252312041e-20
Maximum of d(z^n*J(n)(z))/dz-z^n*J(n-1)(z) =  4.45714481997657e-24
Maximum of d(z^(-m)*J(m)(z))/dz+z^(-m)*J(m+1)(z) =  3.682714599394252e-33
Maximum of z*dJ(m)(z)/dz-m*J(m)+z*J(m+1) =  1.4907779871675686e-19
Enter the degree of the polynomial:9
Enter the total number of roots:10
---------------------------------------------------------
 a=			  b=		  	I=
---------------------------------------------------------
17.241220		 13.354300		 -0.000042
17.241220		 17.241220		 1.000051
17.241220		 20.807048		 -0.000057
17.241220		 24.233885		 0.000062
17.241220		 27.583749		 -0.000067
17.241220		 30.885379		 0.000071
17.241220		 34.154378		 -0.000075
17.241220		 37.400100		 0.000078
17.241220		 40.628554		 -0.000081
17.241220		 43.843801		 0.000084
