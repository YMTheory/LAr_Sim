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
    lr = 0
    lL = 0

    def __init__(self, m_lr, m_lL):
        self.lr = m_lr
        self.lL = m_lL
        LArRindex.LoadData()
        LArTrans.LoadData()

    def __call__(self, *par):
        a0      = par[0]
        aUV     = par[1]
        aIR     = par[2]
        T_v     = par[3]     # temperature of group velocity measurements
        T_t     = par[4]     # temperature of transmission experiment
        rhop0   = par[5]
        rhop1   = par[6]
        delta   = par[7]
        k0      = par[8]
        k1      = par[9]
        R       = par[10]
        mu1     = par[11]
        sigma1  = par[12]
        A2      = par[13]
        mu2     = par[14]
        sigma2  = par[15]
        nu_f    = par[16]

        LArRindex.seta0(a0)
        LArRindex.setaUV(aUV)
        LArRindex.setaIR(aIR)
        LArRindex.setrhop0(rhop0)
        LArRindex.setrhop1(rhop1)

        LArGroupVelocity.setT(T_v)

        LArTrans.setdelta(delta)
        LArTrans.setT(T_t)
        LArTrans.setk0(k0)
        LArTrans.setk1(k1)
        LArTrans.setR(R) 
        LArTrans.setmu1(mu1)
        LArTrans.setsigma1(sigma1)
        LArTrans.setA2(A2)
        LArTrans.setmu2(mu2)
        LArTrans.setsigma2(sigma2)
        LArTrans.setnuf(nu_f)

        chi2 = 0
        chi2 += LArRindex.GetChi2_new()
        chi2 += LArRindex.GetPulls()

        chi2 += LArGroupVelocity.GetChi2()
        chi2 += LArGroupVelocity.GetPulls()
        chi2 += self.lr * LArRindex.rindex_func(0.128)    #### Lagrange multiplier

        chi2 += LArTrans.GetChi2()
        chi2 += LArTrans.GetPulls()
        #chi2 += self.lL * LArTrans.lray_func(0.128)    #### Lagrange multiplier


        #print(LArRindex.GetChi2(), LArGroupVelocity.GetChi2(), LArTrans.GetChi2(), chi2)
        return chi2



class LArFitter(object):

    lr = 0
    lL = 0
    verboseLevel = 0
    chi2min = 0

    @staticmethod
    def setlr(val):
        LArFitter.lr = val

    @staticmethod
    def setlL(val):
        LArFitter.lL = val

    @staticmethod
    def setverbose(val):
        LArFitter.verboseLevel = val

    @staticmethod
    def getchi2min():
        return LArFitter.chi2min


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
    def transmission_pdf(x, a0, aUV, aIR, T, rho0, rho1, delta, k0, k1, R, mu1, sigma1, A2, mu2, sigma2):

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
        LArTrans.setR(R)
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
        m = Minuit(least_squares, a0=0.3347, aUV=0.0994, aIR=0.008, delta=0.307, R=0.937, mu1=0.126, sigma1=0.001, A2=0.4, mu2=0.140, sigma2=0.00154, T=86, rho0=-1.6e-4, rho1=0.0487, k0=6.07e-11, k1=-3.17e-9)
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
    def draw_profile1d(x, y, par):
        fig ,ax = plt.subplots()
        ax.plot(x, y, "-", lw=2, color="coral")
        ax.set_xlabel(par, fontsize=14)
        ax.set_ylabel(r"$\Delta\chi^2$", fontsize=14)
        plt.tight_layout()
        plt.savefig("profile1d_%s.pdf"%par)

    @staticmethod
    def draw_cov(mat):
        fig, ax = plt.subplots()
        ax.matshow(mat)
        plt.savefig("covariance_matrix.pdf")


    @staticmethod
    def fit_generic():
        lsq = MyLeastSquares(LArFitter.lr, LArFitter.lL)
        lsq.func_code = make_func_code(['a0', 'aUV', 'aIR', 'T_v', 'T_t', 'rhop0', 'rhop1', 'delta', 'k0', 'k1', 'R', 'mu1', 'sigma1', 'A2', 'mu2', 'sigma2', 'nuf'])
        m = Minuit(lsq, a0=0.3347, aUV=0.0994, aIR=0.008, delta=0.307, R=0.937, mu1=0.127, sigma1=0.001, A2=0.4, mu2=0.140, sigma2=0.00154, T_v=89, T_t=87, rhop0=-1.6e-4, rhop1=0.0487, k0=6.07e-11, k1=-3.17e-9, nuf=0)
        m.errordef=Minuit.LEAST_SQUARES

        m.migrad()
        m.hesse()
        LArFitter.chi2min = m.fval

        if LArFitter.verboseLevel > 0:
            print("===== Fitting Results =====")
            print(m.fmin)
            print(m.values)
            print(m.errors)
            print(m.covariance)
            print("")
            print("===========================================")
            print("")
            LArRindex.setT(LArGroupVelocity.getT())
            print("Refractive index of LAr at 128nm : %.3f" %LArRindex.rindex_func(0.128))
            print("Group velocity of LAr at 128nm : %.3f with measured value %.3f" %(LArGroupVelocity.Calculate(), LArGroupVelocity.getDataY()) )
            print("Rayleigh scattering lenght of LAr at 128nm : %.3f cm" %LArTrans.lray_func(0.128))
            print("Fitting depolarization ratio value: %.3f" %LArTrans.getdelta())
            print("Fitting group velocity temperature: %.2f K"%LArGroupVelocity.getT())
            print("Fitting transmission temperature: %.2f K"%LArTrans.getT())
            print("Parameters in refractive index formula: a0=%.3f, aUV=%.3f, aIR=%.3f"%(LArRindex.geta0(), LArRindex.getaUV(), LArRindex.getaIR()))
            print("Final chi2 = %.2f (Refractive index) + %.2f (group velocity) + %.2f (transmission)"%(LArRindex.GetChi2_new(), LArGroupVelocity.GetChi2(), LArTrans.GetChi2()))
            print("")
            print("===========================================")


            LArRindex.Plot()
            LArTrans.Plot()
        
            ### Draw 1D profile:
            par = "delta"
            mnp = m.mnprofile("delta", bound=(0, 0.5), subtract_min=True)
            LArFitter.draw_profile1d(mnp[0], mnp[1], par)
    
            cov = m.covariance
            LArFitter.draw_cov(cov)


















