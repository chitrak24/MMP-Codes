"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Solving Laplace's Equation
Author: Chitrak Roychowdhury
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#setting maximum iteration, system dimension and grid
max_iter = 500
X = 50
Y = 50
dx = 1

#Boundary Conditions
Pleft = 100
Pright = 75
Ptop = 250
Pbottom = 500

#setting meshgrid
n1, n2 = (50,50)
dimX = np.linspace(0,50,n1)
dimY = np.linspace(0,50,n2)
Xmesh, Ymesh = np.meshgrid(dimX,dimY)

P = np.zeros((X,Y))
T  = np.empty((X,Y))

P[(Y-1):,:     ] = Ptop
P[:1    ,:     ] = Pbottom
P[:     ,(X-1):] = Pright
P[:     ,:1    ] = Pleft

# Start the iteration
for iteration in range(0, max_iter):
   for i in range(1, X-1, dx):
       for j in range(1, Y-1, dx):
           P[i,j] = 0.25*(P[i+1][j] + P[i-1][j] + P[i][j+1] + P[i][j-1])
           
print("Iteration Completed")

#Ploting
plt.title("Temperature Cantour")
plt.contourf(Xmesh, Ymesh, P, 500, cmap=cm.seismic) # interpolating points = 100
plt.text(21, 4,r'$\phi_{bot} = $'  +str(Pbottom), size = 20)
plt.text(21,45,r'$\phi_{top} = $'  +str(Ptop), size = 20)
plt.text(1, 25,r'$\phi_{left} = $' +str(Pleft),size = 20)
plt.text(42,25,r'$\phi_{right} = $'+str(Pright),size = 20)
plt.colorbar()
plt.grid()
plt.savefig('laplace2D_sol.pdf')
plt.show()
