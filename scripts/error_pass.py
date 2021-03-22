# calculate error pass in rindex and LRay

import numpy as np

def rindex_calc(fitP):
    a0, aUV, aIR, delta = fitP
    l, lUV, lIR = 128, 106.6, 908.3  # nm
    A = a0 + aUV*l**2/(l**2-lUV**2) + aIR*l**2/(l**2-lIR**2)
    n = np.sqrt(1+3*A/(3-A))
    return n

def rindex_error(fitP, fitE, fitC):
    # elements in fitC: c12 ,c13, c23
    a0, aUV, aIR, delta = fitP
    a0e, aUVe, aIRe = fitE
    l, lUV, lIR = 128, 106.6, 908.3  # nm
    p0 = 1
    p1 = l**2/(l**2-aUV**2)
    p2 = l**2/(l**2-aIR**2)
    c01, c02, c12 = fitC
    sigma2 = (a0e*p0)**2 + (aUVe*p1)**2 + (aIRe*p2)**2 \
           + 2*c01*a0e*aUVe*p0*p1 \
           + 2*c02*a0e*aIRe*p0*p2 \
           + 2*c12*a0e*aUVe*p1*p2

    return np.sqrt(sigma2)

def rayleigh_calc(fitP):
    a0, aUV, aIR, delta = fitP
    l, lUV, lIR = 0.128, 0.1066, 0.9083  # nm
    A = a0 + aUV*l**2/(l**2-lUV**2) + aIR*l**2/(l**2-lIR**2)
    n = np.sqrt(1+3*A/(3-A))

    kT = 2.24442E-9
    kB = 1.380649E-23
    T = 90
    f = 1e22

    rayL = 1 / ( 8*np.pi**3/3/l**4 * ((n**2-1)*(n**2+2)/3)**2 * kT * kB * T * f * (6+3*delta)/(6-7*delta) )    

    return rayL


def main():
    fitP = [0.335, 0.099, 0.008, 2.39123e-01]
    rindex_128nm = rindex_calc(fitP)
    print("========> rindex @128nm is %.3f" %rindex_128nm )

    fitE = [0.003, 0.003, 0.003]
    fitC = [-0.614, -0.773, 0.160]
    rindexErr_128nm = rindex_error(fitP, fitE, fitC)
    print("========> rindexErr @128nm is %.4f" %rindexErr_128nm)

    rayL_128nm = rayleigh_calc(fitP)
    print("========> rayL 128nm is %.4f" %rayL_128nm)

if __name__=='__main__':
    main()