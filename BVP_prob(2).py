"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008
Description  : BVP d^2y/dt^2 = -g ; y(0)=0; y(10)=50; (Rocket Problem) 
Author       : Chitrak Roychowdhury
"""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# Feed system size and fixed step-size 
L = 10; dt = 5.0/L

# Construct finite (central) difference A matrix and b
t = np.linspace(0, 5, L+1)
print(t)
b = np.zeros((L+1, 1)).ravel()
print(b)
b[1:L] = -9.8 * pow(dt,2)
b[-1]=50
print(b)
main_diag = -2*np.ones((L+1,1)).ravel()
print("main-diag",main_diag)
off_diag =   1*np.ones((L,  1)).ravel()
print("off-diag",off_diag)
a = main_diag.shape[0]
print(a)
diagonals = [main_diag, off_diag, off_diag]
print("\n diagonals ",diagonals)
A = sparse.diags(diagonals, [0,-1,1]).toarray()

# Enforce y(1) = 0; y(L+1) = 50
A[0, 0  ] = 1
A[0, 1  ] = 0
A[L, L  ] = 1
A[L, L-1] = 0

print(A)
# Solve y = A^-1 * b
y = np.linalg.solve(A,b)

# To match with exact solution
tf = np.linspace(0, 5, 10*L)
yexact = -4.9*pow(tf,2) + 34.5*tf

# Ploting
plt.figure()
plt.plot(t,  y,      'D', c='y',   label= 'FD approximation',               lw=3, ms=8)
plt.plot(tf, yexact, '-', c='r',   label=r'Exact Solution $y=34.5t-4.9*t^2$',        lw=3, ms=8)
plt.legend(loc='best', prop={'size':18})
plt.title(r'BVP using FD Method : $\frac{d^2y}{dt^2} = -g$',size=18)
plt.xlabel('t', size = 16); plt.xticks(size = 14);
plt.ylabel('y', size = 16); plt.yticks(size = 14)
plt.grid();
plt.savefig('bvp(2).pdf') 
plt.show()
