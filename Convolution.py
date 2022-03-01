"""
Registration: 012-1111-0461-20;
Roll: 203012-21-0008
Description: Convolution of Two Functions
Author: Chitrak Roychowdhury
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

conv1 =1; conv2 =1; conv3=1;

#~~~~~~~~~~Convolution of two Gaussian Functions~~~~~~~~~~

if(conv1):
    def G(x,mu,sig):
        return np.exp(-((x-mu)**2)/(2*sig**2))/np.sqrt(2*np.pi*sig**2)

    #Chosing two normal distributions for convolution

    mu1 = 0; sig1 = 0.3;
    mu2 = 0; sig2 = 0.6;
    x = np.linspace(-2,2,1000)
    dx = x[1] - x[0]

    #Performing the convolution
    convolution = np.convolve(G(x,mu1,sig1), G(x,mu2,sig2), mode='same')*dx

    sigc = np.sqrt(sig1**2 + sig2**2)              #convolution std
    ampc = 1/np.sqrt(2*np.pi*(sig1**2 + sig2**2))  #convolution amplitude

    #plot
    plt.figure()
    plt.plot(x, G(x,mu1,sig1), '--b+', lw='4', label= r"$\mathcal{N}_{1}("+str(mu1)+","+str(sig1)+")$")
    plt.plot(x, G(x,mu2,sig2), '--g<', lw='4', label= r"$\mathcal{N}_{2}("+str(mu2)+","+str(sig2)+")$")
    plt.plot(x, convolution, '--rx', lw='4')
    plt.title("(Convolution) Amplitude = "+str(round(ampc,2)))
    plt.legend(loc='best')
    plt.xlabel("$x \longrightarrow $", size='16')
    plt.xticks()
    plt.ylabel(r"$\mathcal{N}(\mu, \sigma)$", size='16')
    plt.yticks()
    plt.grid(True)
    plt.axhline(0, c='brown', ls=':', lw='2')
    plt.axvline(0, c='brown', ls=':', lw='2')
    plt.savefig("conv_gauss.pdf")
    plt.show()

#~~~~~~~~~~Convolution of e^(-x) and sin(x)~~~~~~~~~~

if(conv2):
    def f(x):
        return np.exp(-x)
    def g(x):
        return np.sin(x)

    x = np.linspace(0,30,1000)
    dx = x[1] - x[0]

    #Performing the convolution
    convolution = np.convolve(f(x), g(x), mode='same')*dx

    #Plot
    plt.figure()
    plt.plot(x, f(x), '-', color='blue', lw='2', label='$e^{-x}$')
    plt.plot(x, g(x), '--', color='red', lw='2', label='$\sin(x)$')
    plt.plot(x, convolution, 'o', color='yellow', lw='1',label='exact')
    plt.title("Convolution of $e^{-x}$ and $\sin(x)$")
    plt.legend(loc='best')
    plt.xlabel("$x \longrightarrow $", size='16')
    plt.xticks()
    plt.yticks()
    plt.grid(True)
    plt.axhline(0, c='brown', ls=':', lw='2')
    plt.axvline(0, c='brown', ls=':', lw='2')
    plt.savefig("conv_e^-x-Sin(x).pdf")
    plt.show()

#~~~~~~~~~~Convolution of two unit Step Functions~~~~~~~~~~

if(conv3):
    x = np.linspace(-10,10,1000)
    dx = x[1] - x[0]

    def f(x):
        if x>0:
            return 1
        else:
            return 0
    f = np.vectorize(f)
    
    #Performing the convolution
    conv=[]
    for i in x:
        p = np.linspace(0,i,100)
        h = f(p)*f(i-p)
        I = simps(h,p)
        conv.append(I)

    #Exact result
    def exact(x):
        if x>0:
            return x
        else:
            return 0
    exact = np.vectorize(exact)
        
    #Plot
    plt.figure()
    plt.plot(x, f(x), '-', color='blue', lw='2', label='Unit Step Func')
    plt.plot(x, exact(x), 'o', color='red', lw='1', label='Convolution Result')
    plt.title("Convolution of Step Functions")
    plt.legend(loc='best')
    plt.xlabel("$x \longrightarrow $", size='16')
    plt.xticks()
    plt.yticks()
    plt.grid(True)
    plt.axhline(0, c='brown', ls=':', lw='2')
    plt.axvline(0, c='brown', ls=':', lw='2')
    plt.savefig("conv_step.pdf")
    plt.show()





    
