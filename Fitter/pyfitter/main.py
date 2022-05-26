import numpy as np
import matplotlib.pyplot as plt

from LArRindex import LArRindex
from LArGroupVelocity import LArGroupVelocity
from LArTrans import LArTrans
from LArFitter import LArFitter

def pretest():
    # initial values from previous fitter
    a0=0.3347 
    aUV=0.0994 
    aIR=0.008 
    delta=0.307 
    A1=0.4 
    mu1=0.126 
    sigma1=0.001 
    A2=0.4 
    mu2=0.14012 
    sigma2=0.00154 
    T=86 
    rho0=-1.6e-4 
    rho1=0.0487 
    k0=6.07e-11 
    k1=-3.17e-9

    LArRindex.LoadData()
    LArRindex.seta0(a0)
    LArRindex.setaUV(aUV)
    LArRindex.setaIR(aIR)
    LArRindex.setT(T)
    LArRindex.setrho0(rho0)
    LArRindex.setrho1(rho1)
    LArRindex.Calculate()
    LArRindex.Plot()
    print(LArRindex.GetChi2())

    print(LArGroupVelocity.getDataY(), LArGroupVelocity.Calculate())

    LArTrans.setdelta(delta)
    LArTrans.setT(T)
    LArTrans.setk0(k0)
    LArTrans.setk1(k1)
    LArTrans.setA1(A1)
    LArTrans.setmu1(mu1)
    LArTrans.setsigma1(sigma1)
    LArTrans.setA2(A2)
    LArTrans.setmu2(mu2)
    LArTrans.setsigma2(sigma2)
    
    LArTrans.LoadData()
    LArTrans.Calculate()
    LArTrans.Plot()
    print(LArTrans.lray_func(0.128), LArTrans.Tabs_func(0.128), LArTrans.fresnel_func(0.128), LArTrans.transmission_func(0.128))
    print(LArTrans.GetChi2())


def check_fresnel():
    dx = np.arange(125, 150, 1)
    dy = LArTrans.fresnel_func(dx/1000.)

    fig, ax = plt.subplots()
    ax.plot(dx, dy, "-", color="coral", lw=2)
    ax.set_xlabel("wavelength [nm]", fontsize=14)
    ax.set_ylabel("Fresnel correction", fontsize=14)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__" :
   
    #pretest()
    #check_fresnel()
    
    LArFitter.initialize()
    #LArFitter.fit()
    LArFitter.fit_generic()


