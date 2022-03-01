"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008
Description  : Improper Integrals
Author       : Chitrak Roychowdhury
"""

import numpy as np
import scipy.integrate as sci
import matplotlib.pyplot as plt
from scipy.signal import gaussian

# Logical case switch for different problems to choose from 
prob1=1; prob2=1; 

#Improper Integral
if(prob1):
    a = 3
    def f(x,mu,sig):
        return np.exp(-(x-mu)**2/(2.0*sig**2))*(x+a)/np.sqrt(2.0*np.pi*sig**2)
    
    # Enter standard deviation and integration limits
    mu = 0.0; sig = 0.3; lowlimit = -np.inf; uplimit = np.inf;

    # Eradicating the Crossover problem by narrowing the bracket
    lowlimit = mu - 10*sig; uplimit = mu + 10*sig;

    # Peform the infinite interval integral using Quadrature
    I, err = sci.quad(f,lowlimit,uplimit,args=(mu,sig))

    # Print result
    print ('Integral computed value = ', I, ' with error = ', err)
    
    #Plot
    X=np.linspace(-10,10,1000)
    plt.figure(1)
    plt.plot(X, f(X,mu,sig),lw='1', marker='o',  color='r', label=r'$I$') 
    plt.legend(loc='best')
    plt.grid(True)
    plt.axis([-10, 10, -10, 10])
    plt.title('Improper Integral', fontsize = 16)
    plt.xlabel('$x \longrightarrow $', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.ylabel(r'$\frac{x+3}{\sqrt{2\pi\sigma^{2}}}e^{\frac{-(x-\mu)^{2}}{2\sigma^{2}}}$',fontsize = 16)
    plt.yticks(fontsize = 14)
    plt.savefig("impintegral.pdf") 
    plt.show()

if(prob2):
    def f(x,a,b,c):
        return np.exp(-a*x**2 + b*x + c)

    #Enter the coefficients and integral limits
    a=1.0;
    b=0.0;
    c=1.0;
    llim = -np.inf;
    ulim = np.inf;

    #Performing integration using Quadrature
    I_num, err= sci.quad(f, llim, ulim, args=(a,b,c))
    I_theo = np.sqrt(np.pi/a)*np.exp(b**2/(4*a) + c)

    #Print
    print ('Integral_',llim,'^',ulim,'e^(-',a,'x^2+',b,'+x','+',c,') dx=' , I_num)
    print ('Theoretical value of Integral :', I_theo)
    print ('Absolute Error', err)
    print ('Relative Error', I_num-I_theo)

    #Plot
    X=np.linspace(-10,10,1000)
    plt.figure(1)
    plt.plot(X, f(X,a,b,c), ls='-', lw='5', color='b', label=r'$G$')
    plt.legend(loc='best')
    plt.grid(True)
    plt.axis([-10, 10, -10, 10])
    plt.title('Gaussian Integral', fontsize = 16)
    plt.xlabel('$x \longrightarrow $', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.ylabel(r'$e^{-ax^{2}+bx+c}$',fontsize = 16)
    plt.yticks(fontsize = 14)
    plt.savefig("gaussintegral.pdf") 
    plt.show()
    

    

 
