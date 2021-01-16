#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

def rindex_ba(wl):
    l = wl
    lUV = 0.1066
    lIR = 0.9083
    a0 = 0.335
    aUV = 0.0987
    aIR = 0.008
    A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ;
    n = np.sqrt(1+3*A/(3-A))
    return n

def rindex_ba_mid(wl):
    l = wl
    lUV = 0.1066
    lIR = 0.9083
    a0 = 0.335
    aUV = 0.0987
    aIR = 0.008
    A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ;
    return A

def rindex_our_mid(l, a, b, c):
    A = 1.2055e-2*2/3.;
    rho_ratio = 34.49/(44.66e-3);
    l1 = 91.012;
    l2 = 89.892;
    l3 = 214.02;
    return A*2/3.*rho_ratio * (a/(l1-1/l/l) + b/(l2-1/l/l) + c/(l3-1/l/l));

def rindex_our(wl, a, b, c):
    A = 1.2055e-2*2/3.;
    rho_ratio = 34.49/(44.66e-3);
    l1 = 91.012;
    l2 = 89.892;
    l3 = 214.02;
    l = wl
    return np.sqrt((3/(1-(A*rho_ratio*(a/(l1-1/l/l)+b/(l2-1/l/l)+c/(l3-1/l/l)))))-2);

def raylength(wl, delta, n):
    l = wl
    kT = 2.24442E-9
    kB = 1.380649E-23
    T = 90
    f = 1e22;
    rayL = 1 / (8*np.pi**3/3/l**4 * ((n**2-1)*(n**2+2)/3)**2 * kT * kB * T * f *(6+3*delta)/(6-7*delta))

    return rayL


def main():
    xarr  = np.arange(0.12, 0.7, 0.001)
    yarr1, yarr2, yarr3 = [], [], []
    for i in xarr:
        yarr1.append(rindex_ba_mid(i))
        yarr2.append(rindex_our_mid(i, 0.2075, 0.0415, 4.3330)*3)
    
    plt.plot(xarr, yarr1, "-" , label="model1: x" )
    plt.plot(xarr, yarr2, "--", label="model2: x" )
    #plt.plot(xarr, yarr3, "-.", label="f0=%.2f,f1=%.2f,f2=%.2f" %(0, 0.338, 4.068) )

    plt.xlabel("wavelength/um"); plt.ylabel(r"$4\pi N \alpha$")
    plt.legend()
    #plt.savefig("compare_100+nm.pdf")
    plt.show()

if __name__ == '__main__':
    plt.style.use("seaborn-deep")

    main()