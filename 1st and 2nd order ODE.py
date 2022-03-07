"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008
Description  : Solution of first & second Order Differential Equations using Odeint
Author       : Chitrak Roychowdhury
"""

import numpy as np
from scipy.integrate import quad
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Logical case switch for different problems to choose from 
gquad=1; radiodec=1; dshm=1; fshm=1; vanderpol=1; lorentz=1;

if(gquad): # Gaussian_Quadrature
    s = lambda x: x**2
    print ('Integral 0 to 2 x^2dx using Gauss Quadrature : ', quad(s,0,2))

if(radiodec):
    print('~~~ 1ST ORDER LINEAR ODE : Radioactive Decay of Nuclear Mass ~~~')
    k = 1.5                   # parameter(Force Constant)
   
    def f(x,t):
        dxdt = -k*x
        return dxdt

    t = np.linspace(0,10,100) # Creating time interval; 100 values in [0-100]
    x0 = 100                  # initial value 
    sol = odeint(f, x0, t)    # solution using odeint

    plt.figure()
    plt.plot(t, sol, 'o', color='r', label='x(t)', lw=2, ms=8)
    plt.plot(t, x0*np.exp(-k*t), 'k-', label=r'$x_0 e^{-kt}$', lw=2, ms=8)
    plt.legend(loc='best')
    plt.grid()
    plt.axis([0, 10, 0, 100])
    plt.title(r'Nuclear Decay Curve $(k=1.5,x_0=100)$', fontsize=16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('$ Time \longrightarrow $', fontsize = 16)
    plt.ylabel('$ Nuclear Mass X(t) $',  fontsize = 16)
    plt.savefig('radiodecay.pdf')
    plt.show()
    
if(dshm):
    print("~~~2ND ORDER LINEAR ODE : Damped SHM d2x/dt2 + lambda*dx/dt + kx = 0~~~")

    k=1 ; lam=0.1  #parameters for Damped Harmonic Oscillation

    def dshm(u,t):
        x = u[0] ; y = u[1]
        dxdt = y
        dydt = -k*x-lam*y
        return np.array([dxdt,dydt])

    u0 = [1,0]     #initial values
    t = np.linspace(0,100,1000)
    sol = odeint(dshm, u0, t)
    x1 = sol[:,0]; y1 = sol[:,1]

    #plot
    plt.figure(2)
    plt.plot(t, x1, '^',  label='x(t)', color='g', lw=2, ms=3)
    plt.plot(t, y1, 'v',  label='v(t)', color='b', lw=2, ms=3)
    plt.legend(loc='best') 
    plt.axis([0, 50, -1, 1])
    plt.grid(True)
    plt.axhline(lw=2) # draw a horizontal line
    plt.suptitle('Damped Harmonic Motion',fontsize=16)
    plt.text(25,-0.3,r'$k=1,\lambda=0.1 $', fontsize=20)
    plt.text(13,-0.65,r'$\frac{d^2x}{dt^2}+\lambda\frac{dx}{dt}+kx=0$', fontsize=20)
    plt.xlabel('$ Time \longrightarrow $', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.ylabel('Displacement', fontsize = 16)
    plt.yticks(fontsize = 14)
    plt.savefig('dampedshm.pdf')
    plt.show()

if(fshm):
    print("~~~2ND ORDER LINEAR ODE : Forced SHM d2x/dt2 + lambda*dx/dt + kx = acos(wt)~~~")

    k=1 ; lam=0.2 ; w=0.3 ; a=0.3 #parameters for Damped Harmonic Oscillation

    def dshm(u,t):
        x = u[0] ; y = u[1]
        dxdt = y
        dydt = -k*x-lam*y-a*np.cos(w*t)
        return np.array([dxdt,dydt])

    u0 = [1,0]     #initial values
    t = np.linspace(0,100,1000)
    sol = odeint(dshm, u0, t)
    x2 = sol[:,0]; y2 = sol[:,1]

    #plot
    plt.figure(2)
    plt.plot(t, x2, 'k-',  label='x(t)' ,  lw=2, ms=3)
    plt.plot(t, y2, 'r.-',  label='v(t)', lw=2, ms=3)
    plt.legend(loc='best') 
    plt.axis([0, 50, -1, 1])
    plt.grid(True)
    plt.axhline(lw=2) # draw a horizontal line
    plt.suptitle('Forced Harmonic Motion',fontsize=16)
    plt.text(25,0.5,r'$k=1,\lambda=0.1,\omega=0.2,a=0.2$', fontsize=15)
    plt.text(13,-0.65,r'$\frac{d^2x}{dt^2}+\lambda\frac{dx}{dt}+kx=acos(\omega t)$', fontsize=15)
    plt.xlabel('$ Time \longrightarrow $', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.ylabel('Displacement', fontsize = 16)
    plt.yticks(fontsize = 14)
    plt.savefig('forcedshm.pdf')
    plt.show()

if(vanderpol):
   print ('~~~ 2ND ORDER NONLINEAR ODE (Vanderpol Oscillator) : d2x/dt2 - mu(1-x^2)*dx/dt + beta*x = 0 ~~~ ')
   print ('~~~ dydt = x/mu & dxdt = mu(x-x^3/3-beta*y); [mu & beta > 0] ~~~ ')


   mu, beta = 1.0, 1.0   # Parameters of nonlinearity & Hookean Elasticity

   def vanderpol(X, t):
       x = X[0]
       y = X[1]
       dxdt = mu*(x - x**3/3.0 - beta*y)
       dydt = x/mu
       return [dxdt, dydt]

   # main
   x0 = [1, 2]     # initial values
   t = np.linspace(0, 9000, 450)
   sol = odeint(vanderpol, x0, t)
   x = sol[:,0]; y = sol[:,1]

   # plot
   plt.figure()
   plt.subplot(2,1,1); # dynamics
   plt.plot(t, x, 'k-',  label='x(t)', lw=2, ms=6)
   plt.plot(t, y, 'r.-', label='y(t)', lw=2, ms=6)
   plt.legend(loc='best') 
   plt.grid()
   plt.suptitle(r'Vanderpol Oscillator: $\frac{d^2x}{dt^2}+\mu(1-x^2)\frac{dx}{dt}+\beta x=0$', fontsize=20)
   plt.title(r'$\mu='+str(mu)+', \eta='+str(beta)+'$', fontsize=20)
   plt.xlabel('Time', fontsize = 16);
   plt.xticks(fontsize = 14)
   plt.ylabel('Displacement', fontsize = 16);
   plt.yticks(fontsize = 14)
   
   plt.subplot(2,1,2); # phase portrait
   plt.plot(x, y,       'go-', label='x(t) vs y(t)', lw=1, ms=6)
   plt.plot(x[0], y[0], 'b*', label='Initial Value', lw=2, ms=20)
   plt.legend(loc='best') 
   plt.xlabel('X Displacement', fontsize = 16);
   plt.xticks(fontsize = 14)
   plt.ylabel('Y Displacement', fontsize = 16);
   plt.yticks(fontsize = 14)
   plt.grid()
   plt.savefig('vanderpol.pdf')
   plt.show()

if(lorentz):
   print ('~~~ 2ND ORDER NONLINEAR ODE : Lorentz Attractor ~~~ ')                       
   print ('~~~ dx/dt=sigma*(y-x), dy/dt=x*(rho-z)-y, dz/dt=x*y-beta*z  ~~~')

   sig, rho, beta = 10.0, 28.0, 8.0/3
   def loratr(u, t):
       x,y,z = u[0],u[1],u[2]
       dxdt = sig*(y-x)
       dydt = x*(rho-z)-y
       dzdt = x*y-beta*z
       return [dxdt, dydt, dzdt]
    
   u0 = [0, 1.0, 0]
   t = np.linspace(0,50,50000)
   sol = odeint(loratr, u0, t)
   x, y, z = sol[:,0], sol[:,1], sol[:,2]

   # Plot
   plt.figure(3)
   plt.subplot(2,2,1)
   plt.plot(x, z, 'r-', label='X-Z', lw=2, ms=8)
   plt.legend(loc='best') 
   plt.axis([-20, 20, 0, 50])
   plt.suptitle('Lorentz Attractor', fontsize=18)
   plt.xlabel('X', fontsize = 16)
   plt.xticks(fontsize = 14)
   plt.ylabel('Z', fontsize = 16)
   plt.yticks(fontsize = 14)

   plt.subplot(2,2,2)
   plt.plot(x, y, 'k-', label='X-Y', lw=2, ms=8)
   plt.legend(loc='best') 
   plt.axis([-20, 20, -30, 30])
   plt.xlabel('X', fontsize = 16)
   plt.xticks(fontsize = 14)
   plt.ylabel('Y', fontsize = 16)
   plt.yticks(fontsize = 14)

   plt.subplot(2,2,3)
   plt.plot(y, z, 'g.', label='Y-Z', lw=2, ms=2)
   plt.legend(loc='best') 
   plt.axis([-30, 30, 0, 50])
   plt.xlabel('Y', fontsize = 16)
   plt.xticks(fontsize = 14)
   plt.ylabel('Z', fontsize = 16)
   plt.yticks(fontsize = 14)
   
   plt.text(50, 35, r'$\frac{dx}{dt}=\sigma(y-x)$', fontsize=20)
   plt.text(50, 20, r'$\frac{dy}{dt}= x(\rho-z)-y$', fontsize=20)
   plt.text(50, 5, r'$\frac{dz}{dt}= xy-\beta z$', fontsize=20)   
   plt.savefig('lorentz.pdf')
   plt.show()
