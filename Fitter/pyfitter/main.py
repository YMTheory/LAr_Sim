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
    T_v = 89
    T_t=86 
    rho0=-1.6e-4 
    rho1=0.0487 
    k0=6.07e-11 
    k1=-3.17e-9

    LArRindex.LoadData()
    LArRindex.seta0(a0)
    LArRindex.setaUV(aUV)
    LArRindex.setaIR(aIR)
    #LArRindex.setT(T)
    LArRindex.setrhop0(rho0)
    LArRindex.setrhop1(rho1)
    LArRindex.Calculate_new()
    LArRindex.Plot()
    print(LArRindex.GetChi2_new())

    LArGroupVelocity.setT(T_v)
    print(LArGroupVelocity.getDataY(), LArGroupVelocity.Calculate(), LArGroupVelocity.GetChi2())

    LArTrans.setdelta(delta)
    LArTrans.setT(T_t)
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


if __name__ == "__main__" :
   
    #pretest()

    LArFitter.setverbose(0)
    LArFitter.setToyMC(True)
    LArFitter.initialize()
    delta_arr, R_arr, chi2_arr = [], [], []
    for i in range(300, 400, 1):
        print("Running toyMC dataset %d"%i)
        LArFitter.setseed(i)
        LArFitter.fit_generic()
        delta_arr.append(LArTrans.getdelta())
        R_arr.append(LArTrans.getR())
        chi2_arr.append(LArFitter.getchi2min())

        LArRindex.setT(LArGroupVelocity.getT())
        print(LArTrans.getdelta(), LArTrans.getR(), LArRindex.rindex_func(0.128), LArTrans.lray_func(0.128), LArFitter.getchi2min())

    ## Lagrange multiplier : refractive index and Rayleigh scattering length
    #rindex, lray, chi2min = [], [], []
    #for i in np.arange(-2000, 2010, 10):
    #    print("factor %.2f"%i)
    #    LArFitter.setlr(i)
    #    LArFitter.fit_generic()
    #    
    #    LArRindex.setT(LArGroupVelocity.getT())
    #    rindex.append(LArRindex.rindex_func(0.128))
    #    #LArRindex.setT(LArTrans.getT())
    #    #lray.append(LArTrans.lray_func(0.128))
    #    chi2min.append(LArFitter.getchi2min())
    #    print(rindex[-1], chi2min[-1])
    #    #print(lray[-1], chi2min[-1])

