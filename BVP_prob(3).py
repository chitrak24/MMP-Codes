"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008
Description  : BVP y'' = x^2-2; f(0)=0; f(3)=2; Exactsol y = x^4/12-x^2+(17/12)*2
Author       : Chitrak Roychowdhury
"""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# Feed system size and fixed step-size 
L = 30; dx = 3.0/L

# Construct finite (central) difference A matrix and b
x = np.linspace(0, 3, L+1)
print(x)
b = np.zeros((L+1, 1)).ravel()
print(b)
b[1:L] = pow(dx,2) * (x[1:L]**2-2)
b[-1]= 2
main_diag = -2*np.ones((L+1,1)).ravel()
print("main-diag",main_diag)
off_diag =   1*np.ones((L,  1)).ravel()
print("off-diag",off_diag)
a = main_diag.shape[0]
print(a)
diagonals = [main_diag, off_diag, off_diag]
print("\n diagonals ",diagonals)
A = sparse.diags(diagonals, [0,-1,1]).toarray()

# Enforce f(1) =0 & f(L+1) = 2
A[0, 0  ] = 1
A[0, 1  ] = 0
A[L, L  ] = 1
A[L, L-1] = 0

print(A)
# Solve y = A^-1 * b
f = np.linalg.solve(A,b)

# To match with exact solution
xf = np.linspace(0, 3, 10*L)
fexact = pow(xf,4)/12 - pow(xf,2) +(17/12)*xf

# Ploting
plt.figure()
plt.plot(x,  f,      '^', c='k', label='FD approximation',          lw=3, ms=8)
plt.plot(xf, fexact, '-', c='b',  label=r'Exact Solution $y=\frac{x^4}{12}-x^2+\frac{17}{12}x$',        lw=3, ms=8)
plt.legend(loc='best', prop={'size':18})
plt.title(r'BVP using FD Method : $\frac{d^2y}{dx^2} = x^2-2$',size=18)
plt.xlabel('x', size = 16); plt.xticks(size = 14);
plt.ylabel('y', size = 16); plt.yticks(size = 14)
plt.grid();
plt.savefig('bvp(3).pdf') 
plt.show()
