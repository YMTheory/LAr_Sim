from iminuit import cost, Minuit
from iminuit.cost import LeastSquares
from iminuit.util import describe, make_func_code

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import numba as nb

from LArRindex import LArRindex
from LArTrans import LArTrans
from LArGroupVelocity import LArGroupVelocity

class MyLeastSquares:
    """
    Generic least squares functions with pull terms
    """
    errordef = Minuit.LEAST_SQUARES # for Minuit to compute errors correctly

    def __init__(self):
        LArRindex.LoadData()
        LArTrans.LoadData()

    def __call__(self, *par):
        a0 = par[0]
        aUV = par[1]
        aIR = par[2]
        T = par[3]
        rho0 = par[4]
        rho1 = par[5]
        delta = par[6]
        k0 = par[7]
        k1 = par[8]
        A1 = par[9]
        mu1 = par[10]
        sigma1 = par[11]
        A2 = par[12]
        mu2 = par[13]
        sigma2 = par[14]

        LArRindex.seta0(a0)
        LArRindex.setaUV(aUV)
        LArRindex.setaIR(aIR)
        LArRindex.setT(T)
        LArRindex.setrho0(rho0)
        LArRindex.setrho1(rho1)

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

        chi2 = 0
        chi2 += LArRindex.GetChi2()
        chi2 += LArGroupVelocity.GetChi2()
        chi2 += LArTrans.GetChi2()
        return chi2





class LArFitter(object):

    @staticmethod
    def initialize():
        LArRindex.LoadData()
        LArTrans.LoadData()

    @staticmethod
    def rindex_pdf(x, a0, aUV, aIR, T, rho0, rho1):
        LArRindex.seta0(a0)
        LArRindex.setaUV(aUV)
        LArRindex.setaIR(aIR)
        LArRindex.setT(T)
        LArRindex.setrho0(rho0)
        LArRindex.setrho1(rho1)

        return LArRindex.rindex_func(x)


    @staticmethod
    def groupvelocity_pdf(x, a0, aUV, aIR, T, rho0, rho1):
        LArRindex.seta0(a0)
        LArRindex.setaUV(aUV)
        LArRindex.setaIR(aIR)
        LArRindex.setT(T)
        LArRindex.setrho0(rho0)
        LArRindex.setrho1(rho1)
        return LArGroupVelocity.Calculate()

    
    @staticmethod
    def transmission_pdf(x, a0, aUV, aIR, T, rho0, rho1, delta, k0, k1, A1, mu1, sigma1, A2, mu2, sigma2):

        LArRindex.seta0(a0)
        LArRindex.setaUV(aUV)
        LArRindex.setaIR(aIR)
        LArRindex.setT(T)
        LArRindex.setrho0(rho0)
        LArRindex.setrho1(rho1)

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
        
        return LArTrans.transmission_func(x)


    @staticmethod
    def fit():
        rindex_wavelength = LArRindex.getDataX()
        rindex_data       = LArRindex.getDataY()
        rindex_err        = LArRindex.getDataYerr()

        least_squares1    = LeastSquares(rindex_wavelength, rindex_data, rindex_err, LArFitter.rindex_pdf)


        groupvelocity_wavelength = LArGroupVelocity.getDataX()
        groupvelocity_data       = LArGroupVelocity.getDataY()
        groupvelocity_err        = LArGroupVelocity.getDataYerr()

        least_squares2           = LeastSquares(groupvelocity_wavelength, groupvelocity_data, groupvelocity_err, LArFitter.groupvelocity_pdf)


        trans_wavelength         = LArTrans.getDataX()
        trans_data               = LArTrans.getDataY()
        trans_err                = LArTrans.getDataYerr()

        least_squares3           = LeastSquares(trans_wavelength, trans_data, trans_err, LArFitter.transmission_pdf)

        
        least_squares = least_squares1 + least_squares2 + least_squares3
        m = Minuit(least_squares, a0=0.3347, aUV=0.0994, aIR=0.008, delta=0.307, A1=0.4, mu1=0.126, sigma1=0.001, A2=0.4, mu2=0.140, sigma2=0.00154, T=86, rho0=-1.6e-4, rho1=0.0487, k0=6.07e-11, k1=-3.17e-9)
        m.limits["a0"] = (0.1, 0.5)
        m.limits["aUV"] = (0.05, 0.15)
        m.limits["aIR"] = (0.001, 0.02)
        m.limits["delta"] = (0, 1.0)


        m.migrad()
        m.hesse()
    
        print("===== Fitting Results =====")
        print(m.fmin)
        print(m.values)
        print(m.errors)
        print(m.covariance)

        LArRindex.Plot()
        LArTrans.Plot()
        
        print("")
        print("===========================================")
        print("Refractive index of LAr at 128nm : %.3f" %LArRindex.rindex_func(0.128))
        print("Rayleigh scattering lenght of LAr at 128nm : %.3f cm" %LArTrans.lray_func(0.128))



    @staticmethod
    def fit_generic():
        lsq = MyLeastSquares()
        lsq.func_code = make_func_code(['a0', 'aUV', 'aIR', 'T', 'rho0', 'rho1', 'delta', 'k0', 'k1', 'A1', 'mu1', 'sigma1', 'A2', 'mu2', 'sigma2'])
        # this fails
        try:
            m = Minuit(lsq, a0=0.3347, aUV=0.0994, aIR=0.008, delta=0.307, A1=0.4, mu1=0.126, sigma1=0.001, A2=0.4, mu2=0.140, sigma2=0.00154, T=86, rho0=-1.6e-4, rho1=0.0487, k0=6.07e-11, k1=-3.17e-9)
            m.errordef=Minuit.LEAST_SQUARES

            print("===== Fitting Results =====")
            print(m.fmin)
            print(m.values)
            print(m.errors)
            print(m.covariance)

            LArRindex.Plot()
            LArTrans.Plot()
        
        except:
            import traceback
            traceback.print_exc()

