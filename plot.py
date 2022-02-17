"""
Registration : 012-1111-0461-20
Roll         : 203012-21-0008 
Description  : Ploting Special Functions
Author       : Chitrak Roychowdhury
"""

import numpy as np
from scipy.special import legendre, hermite, jn, laguerre, chebyt
import matplotlib.pyplot as plt

# Logical case switch for different problems to choose from 
legendr=1; herm=1; bessel=1; laug=1; cheby=1;

if(legendr):
     print("~-~-~-~-~-~-LEGENDRE POLYNOMIALS~-~-~-~-~-~-")
     print("P1(x)",legendre(1), "\nP2(x)",legendre(2), "\nP3(x)",legendre(3), "\nP4(x)",legendre(4), "\nP5(x)",legendre(5), "\nP6(X)",legendre(6), "\nP7(x)",legendre(7), "\nP8(x)",legendre(8), "\nP9(x)",legendre(9), "\nP10(x)",legendre(10))
     x = np.arange(-1,1,0.01)
     p1 = legendre(1); p2 = legendre(2); p3 = legendre(3); p4 = legendre(4); 
     p5 = legendre(5); p6 = legendre(6); p7 = legendre(7); p8 = legendre(8);
     p9 = legendre(9); p10 = legendre(10)

     plt.figure(1)
     plt.plot(x, p1(x), lw=3, ls='-',  color='k', label=r'$P_1$')
     plt.plot(x, p2(x), lw=3, ls='--', color='r', label=r'$P_2$')
     plt.plot(x, p3(x), lw=3, ls='-.', color='g', label=r'$P_3$')
     plt.plot(x, p4(x), lw=3, ls=':',  color='m', label=r'$P_4$')
     plt.plot(x, p5(x), lw=3, ls='-',  color='b', label=r'$P_5$')
     plt.plot(x, p6(x), lw=3, ls='dashdot',  color='k', label=r'$P_6$')
     plt.plot(x, p7(x), lw=3, ls=':',  color='g', label=r'$P_7$')
     plt.plot(x, p8(x), lw=3, ls='--',  color='r', label=r'$P_8$')
     plt.plot(x, p9(x), lw=3, ls='dashed', color='b', label=r'$P_9$')
     plt.plot(x, p10(x), lw=3, ls='-', color='m', label=r'$P_{10}$')

     plt.legend(loc='best')
     plt.grid()
     plt.axis([-1, 1, -1, 1])
     plt.title('Legendre Polynomials', fontsize = 16)
     plt.xlabel('$x \longrightarrow $', fontsize = 16)
     plt.xticks(fontsize = 14)
     plt.ylabel(r'$P_n(x) = \frac{1}{2^n n!}\frac{d^n}{dx^n}(x^2-1)^n$',fontsize = 16)
     plt.yticks(fontsize = 14)
     plt.savefig("legendreplot.pdf") 
     plt.show()

if(herm):
    print("~-~-~-~-~-~-HERMITE POLYNOMIALS~-~-~-~-~-~-")
    print("H1(x)",hermite(1), "\nH2(x)",hermite(2), "\nH3(x)",hermite(3), "\nH4(x)",hermite(4), "\nH5(x)",hermite(5))
    x = np.arange(-10,10,0.1)
    h1 = hermite(1); h2 = hermite(2); h3 = hermite(3); h4 = hermite(4); 
    h5 = hermite(5);
    
    plt.figure(2)
    plt.plot(x, h1(x), lw=2, ls=':', color='k', label=r'$H_1$')
    plt.plot(x, h2(x), lw=2, ls='solid', color='r', label=r'$H_2$')
    plt.plot(x, h3(x), lw=2, ls='-.', color='g', label=r'$H_3$')
    plt.plot(x, h4(x), lw=2, ls='-', color='b', label=r'$H_4$')
    plt.plot(x, h5(x), lw=2, ls=':', color='m', label=r'$H_5$')    
    plt.legend(loc='best')
    plt.grid()
    plt.axis([-10, 10, -8000, 8000])
    plt.title('Hermite Polynomials', fontsize = 16)
    plt.xlabel('$x \longrightarrow $', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.ylabel(r'$H_n(x) = (-1)^ne^{x^2}\frac{d^n}{dx^n}e^{-x^2}$',fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.savefig('hermiteplot.pdf') 
    plt.show()
    
if(bessel):
    print("~-~-~-~-~-~-BESSEL FUNCTION~-~-~-~-~-~-")
    x = np.linspace(0,20,500)
    print("J1(x)",jn(1,5), "\nJ2(x)",jn(2,5), "\nJ3(x)",jn(3,5))
    j0 = jn(0,x); j1 = jn(1,x); j2 = jn(2,x); j3 = jn(3,x);

    plt.figure(3)
    plt.plot(x, j0, lw=2,  marker='*', color='k', label=r'$J_0$')
    plt.plot(x, j1, lw=2,  marker='p', color='r', label=r'$J_1$')
    plt.plot(x, j2, lw=2,  marker='s', color='b', label=r'$J_2$')
    plt.plot(x, j3, lw=2,  marker='8',color='g', label=r'$J_3$')
    plt.legend(loc='best')
    plt.grid()
    plt.axis([0, 20, -1, 1])
    plt.title('Bessel Function: $1^{st}$ kind', fontsize = 16)
    plt.xlabel('$x \longrightarrow $', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.ylabel(r'$J_n(x) = \sum_{m=0}^\infty \frac{(-1)^m}{m!\Gamma(m+n+1)}(\frac{x}{2})^{2m+n}$',fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.savefig('besselplot.pdf') 
    plt.show()
    
if(laug):
    print("~-~-~-~-~-~-LAGUERRE FUNCTION~-~-~-~-~-~-")
    x = np.arange(-10,10,0.01)
    print("\nL0(x)",laguerre(0), "L1(x)",laguerre(1), "\nL2(x)",laguerre(2), "\nL3(x)",laguerre(3), "\nL4(x)",laguerre(4), "\nL5(x)",laguerre(5))
    L0 = laguerre(0); L1 = laguerre(1); L2 = laguerre(2); L3 = laguerre(3); L4 = laguerre(4); L5 = laguerre(5)

    plt.figure(3)
    plt.plot(x, L0(x),  marker='+',  color='k', label=r'$L_0$')
    plt.plot(x, L1(x),  marker='D', color='r', label=r'$L_1$')
    plt.plot(x, L2(x),  marker='d', color='b', label=r'$L_2$')
    plt.plot(x, L3(x),  marker='x', color='g', label=r'$L_3$')
    plt.plot(x, L4(x),  marker='h', color='m', label=r'$L_4$')
    plt.plot(x, L5(x),  marker='H', color='k', label=r'$L_5$') 
    plt.legend(loc='best')
    plt.grid()
    plt.axis([-10, 10, -10, 10])
    plt.title('Laguerre Function:', fontsize = 16)
    plt.xlabel('$x \longrightarrow $', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.ylabel(r'$L_n(x) = \sum_{m=0}^{n} \frac{(-1)^{m} n! x^{m}}{(n-m)! m! m!}$',fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.savefig('laguerreplot.pdf') 
    plt.show()
    
if(cheby):
     print("~-~-~-~-~-~-CHEBYSHEV FUNCTION~-~-~-~-~-~-")
     x = np.linspace(-1,1,1000)
     print("\nT0(x)",chebyt(0), "T1(x)",chebyt(1), "\nT2(x)",chebyt(2), "\nT3(x)",chebyt(3), "\nT4(x)",chebyt(4), "\nT5(x)",chebyt(5))
     t0 = chebyt(0); t1 = chebyt(1); t2 = chebyt(2); t3 = chebyt(3); t4 = chebyt(4); t5 = chebyt(5)

     plt.figure(3)
     plt.plot(x, t0(x), lw=2, ls=':',  color='k', label=r'$T_0$')
     plt.plot(x, t1(x), lw=2, ls='--', color='r', label=r'$T_1$')
     plt.plot(x, t2(x), lw=2, ls='-', color='b', label=r'$T_2$')
     plt.plot(x, t3(x), lw=2, ls='-.', color='g', label=r'$T_3$')
     plt.plot(x, t4(x), lw=2, ls='-', color='m', label=r'$T_4$')
     plt.plot(x, t5(x), lw=2, ls=':', color='k', label=r'$T_5$') 
     plt.legend(loc='best')
     plt.grid()
     plt.axis([-1, 1, -1, 1])
     plt.title('Chebyshev Function:', fontsize = 16)
     plt.xlabel('$x \longrightarrow $', fontsize = 16)
     plt.xticks(fontsize = 14)
     plt.ylabel(r'$T_n(x) = \sum_{m=0}^{\infty} \frac{(-2)^{n}n!}{(2n)!}\sqrt{1-x^2}\frac{d^n}{dx^n}(1-x^2)^{n-\frac{1}{2}}$',fontsize = 14)
     plt.yticks(fontsize = 14)
     plt.savefig('chebyshevplot.pdf') 
     plt.show()
    

     
    
