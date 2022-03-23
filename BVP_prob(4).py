"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008
Description  : BVP d^y/dt^2 + w^2*y = 0; y(0)=1 and y'(0)=0
Author       : Chitrak Roychowdhury
"""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# Feed system size and fixed step-size 
L = 50; h = 10.0/L; w = 1;

# Construct finite (central) difference A matrix and b
t = np.linspace(0, 10, L+1)
print(t)
b = np.zeros((L+1, 1)).ravel()
b[0]=1
print(b)
main_diag = -2 * np.ones((L+1,1)) + pow(w,2)*pow(h,2)
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
tf = np.linspace(0, 10, 10*L)
fexact = np.cos(w*tf)

# Ploting
plt.figure()
plt.plot(t,  f,      's', c='r', label='FD approximation',                     lw=3, ms=8)
plt.plot(tf, fexact, '-.', c='b',  label=r'Exact Solution $y=cos(wt)$',        lw=3, ms=8)
plt.legend(loc='best', prop={'size':18})
plt.title(r'BVP using FD Method : $\frac{d^2y}{dt^2}+w^2y = 0 $',size=18)
plt.xlabel('t', size = 16); plt.xticks(size = 14);
plt.ylabel('y', size = 16); plt.yticks(size = 14)

plt.grid();
plt.savefig('bvp(4).pdf') 
plt.show()
