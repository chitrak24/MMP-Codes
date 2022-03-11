"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008
Description  : BVP y'' = 12x^2; y(0)=0; y(1)=0; Exactsol y = x^4-x
Author       : Chitrak Roychowdhury
"""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# Feed system size and fixed step-size 
L = 100; dx = 1.0/L

# Construct finite (central) difference A matrix and b
x = np.linspace(0, 1, L+1)
print(x)
b = np.zeros((L+1, 1)).ravel()
print(b)
b[1:L] = 12 * pow(dx,2) * (x[1:L]**2)

main_diag = -2*np.ones((L+1,1)).ravel()
print("main-diag",main_diag)
off_diag =   1*np.ones((L,  1)).ravel()
print("off-diag",off_diag)
a = main_diag.shape[0]
print(a)
diagonals = [main_diag, off_diag, off_diag]
print("\n diagonals ",diagonals)
A = sparse.diags(diagonals, [0,-1,1]).toarray()

# Enforce y(1) = y(L+1) = 0
A[0, 0  ] = 1
A[0, 1  ] = 0
A[L, L  ] = 1
A[L, L-1] = 0

print(A)
# Solve y = A^-1 * b
y = np.linalg.solve(A,b)

# To match with exact solution
xf = np.linspace(0, 1, 10*L)
yexact = pow(xf,4) - xf

# Ploting
plt.figure()
plt.plot(x,  y,      'o', c='g', label='FD approximation',          lw=3, ms=8)
plt.plot(xf, yexact, '--', c='b',  label=r'Exact Solution $y=x^4-x$',        lw=3, ms=8)
plt.legend(loc='best', prop={'size':18})
plt.title(r'BVP using FD Method : $\frac{d^2y}{dx^2} = 12x^2$',size=18)
plt.xlabel('x', size = 16); plt.xticks(size = 14);
plt.ylabel('y', size = 16); plt.yticks(size = 14)
plt.grid();
plt.savefig('bvp(1).pdf') 
plt.show()
