"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008
Description  : BVP d^y/dx^2 + 4*y - 4*x = 0; y(0)=0 and y'(0)=0
Author       : Chitrak Roychowdhury
"""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# Feed system size and fixed step-size 
L = 50; h = (np.pi/2-0)/L; 

# Construct finite (central) difference A matrix and b
x = np.linspace(0, np.pi/2, L+1)
print(x)
b = np.zeros((L+1, 1)).ravel()
b[1:L+1] = 4 * pow(h,2) *(x[1:L+1])
print(b)
main_diag = -2 * np.ones((L+1,1)) + 4*pow(h,2)
md= main_diag.ravel()
print("main-diag",md)
off_diag =   1*np.ones((L,  1)).ravel()
print("off-diag",off_diag)
diagonals = [md, off_diag, off_diag]
print("\n diagonals ",diagonals)
A = sparse.diags(diagonals, [0,-1,1]).toarray()

# Enforce boundary conditions
A[0, 0  ] = 1
A[0, 1  ] = 0
A[L, L-1] = 2
print(A)

# Solve y = A^-1 * b
f = np.linalg.solve(A,b)

# To match with exact solution
xf = np.linspace(0, np.pi/2, 10*L)
fexact = xf - np.sin(2*xf)

# Ploting
plt.figure()
plt.plot(x,  f,      'd',  c='c',  label='FD approximation',                       lw=3, ms=8)
plt.plot(xf, fexact, '--', c='r',  label=r'Exact Solution $y=x-sin(2x)$',        lw=3, ms=8)
plt.legend(loc='best', prop={'size':18})
plt.title(r'BVP using FD Method : $\frac{d^2y}{dx^2}+4y-4x = 0 $',size=18)
plt.xlabel('t', size = 16); plt.xticks(size = 14);
plt.ylabel('y', size = 16); plt.yticks(size = 14)

plt.grid();
#plt.savefig('bvp(4).pdf') 
plt.show()
