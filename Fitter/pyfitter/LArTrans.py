import numpy as np
from LArRindex import LArRindex

class LArTrans(object):

    wavelength = []
    trans_data = []
    trans_err = []
    trans_calc = []

    delta = 0
    T = 86
    k0 = 6.0712e-11
    k1 = -3.1699e-09

    A1  = 0.4
    mu1 = 126.51
    sigma1 = 1
    A2 = 0.4
    mu2 = 140.12
    sigma2 = 1.537

    chi2 = 0

    @staticmethod
    def LoadData():
        XeDopedFile = "../data/G140ppb.txt"
        with open(XeDopedFile) as f:
            for lines in f.readlines():
                line = lines.strip("\n")
                data = line.split(" ")
                if float(data[0]) > 125:
                    LArTrans.wavelength.append(float(data[0])/1000.)
                    LArTrans.trans_data.append(float(data[1]))
                    LArTrans.trans_err.append(float(data[2]))



    @staticmethod
    def lray_func_iso(l):
        rindex = LArRindex.rindex_func(l)
        delta = 0
        T = LArTrans.T
        p0 = LArTrans.k0
        p1 = LArTrans.k1
        kT = p0 * T + p1

        kB = 1.380649e-23
        f = 1e22

        pi = np.pi

        rayL = 1 / (8*pi**3/3/l**4 * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));
        return rayL


        
    @staticmethod
    def lray_func(l):
        rindex = LArRindex.rindex_func(l)
        delta = LArTrans.delta
        T = LArTrans.T
        p0 = LArTrans.k0
        p1 = LArTrans.k1
        kT = p0 * T + p1

        kB = 1.380649e-23
        f = 1e22

        pi = np.pi

        rayL = 1 / (8*pi**3/3/l**4 * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));
        return rayL


    @staticmethod
    def Tabs_func(l):
        A1     = LArTrans.A1
        mu1    = LArTrans.mu1
        sigma1 = LArTrans.sigma1
        A2     = LArTrans.A2
        mu2    = LArTrans.mu2
        sigma2 = LArTrans.sigma2

        A_abs = A1*np.exp(-(l-mu1)*(l-mu1)/2/sigma1/sigma1) + A2*np.exp(-(l-mu2)*(l-mu2)/2/sigma2/sigma2);
        T_abs = np.exp( -A_abs*np.log(10.) );
        return T_abs

    @staticmethod
    def fresnel_func(l):
        wl = l * 1000
        E = 1240./wl
        a = 38
        E0 = 12.8234
        gamma = 0.42357
        b = 247.188
        E1 = 18.8448
        
        n_MgF2 = np.sqrt(1+a*(E0*E0-E*E)/((E0*E0-E*E)*(E0*E0-E*E)+gamma*E*E) + b/(E1*E1-E*E))

        n_vac = 1.
        n_LAr = LArRindex.rindex_func(l)

        f_F = 1 / (1-((n_vac - n_MgF2)/(n_vac + n_MgF2))**2/1-((n_MgF2 - n_LAr)/(n_MgF2+n_LAr))**2)**2
        return f_F


    @staticmethod
    def transmission_func(l):
        d = 5.8
        rayL = LArTrans.lray_func(l)
        T_Ray = np.exp(-d/rayL)
        f_F = LArTrans.fresnel_func(l)
        T_abs = LArTrans.Tabs_func(l)

        trans_pred = T_Ray * T_abs * f_F
        return trans_pred


    @staticmethod
    def Calculate():
        LArTrans.trans_calc = []
        for i in LArTrans.wavelength:
            trans_pred = LArTrans.transmission_func(i)
            LArTrans.trans_calc.append(trans_pred)


    @staticmethod
    def GetChi2():
        LArTrans.Calculate()
        for i in range(len(LArTrans.wavelength)):
            LArTrans.chi2 += (LArTrans.trans_data[i] - LArTrans.trans_calc[i])**2 / (LArTrans.trans_err[i])**2
        return LArTrans.chi2

    @staticmethod
    def Plot():
        import matplotlib.pyplot as plt
        LArTrans.Calculate()
        fig, ax = plt.subplots()
        ax.errorbar(LArTrans.wavelength, LArTrans.trans_data, yerr=LArTrans.trans_err, fmt="o", label="Simulation")
        ax.plot(LArTrans.wavelength, LArTrans.trans_calc, "o", label="Calculation")
        ax.legend(prop={"size":14})
        ax.set_xlabel("wavelength [um]", fontsize=14)
        ax.set_ylabel("transmission", fontsize=14)
        ax.grid(True)
        plt.tight_layout()
        plt.savefig("trans.pdf")
                


    ########## Getter Functions ##########

    @staticmethod
    def getdelta():
        return LArTrans.delta

    @staticmethod 
    def getT():
        return LArTrans.T

    @staticmethod
    def getk0():
        return LArTrans.k0


    @staticmethod
    def getk1():
        return LArTrans.k1

    @staticmethod
    def getA1():
        return LArTrans.A1

    @staticmethod
    def getmu1():
        return LArTrans.mu1

    @staticmethod
    def getsigma1():
        return LArTrans.sigma1

    @staticmethod
    def getA2():
        return LArTrans.A2

    @staticmethod
    def getmu2():
        return LArTrans.mu2

    @staticmethod
    def getsigma2():
        return LArTrans.sigma2

    @staticmethod
    def getDataX():
        return LArTrans.wavelength

    @staticmethod
    def getDataY():
        return LArTrans.trans_data

    @staticmethod
    def getDataYerr():
        return LArTrans.trans_err


    ######### Setter Functions ##########
    @staticmethod
    def setdelta(val):
        LArTrans.delta = val

    @staticmethod
    def setT(val):
        LArTrans.T = val

    @staticmethod
    def setk0(val):
        LArTrans.k0 = val

    @staticmethod
    def setk1(val):
        LArTrans.k1 = val

    @staticmethod
    def setA1(val):
        LArTrans.A1 = val

    @staticmethod
    def setmu1(val):
        LArTrans.mu1 = val

    @staticmethod
    def setsigma1(val):
        LArTrans.sigma1 = val

    @staticmethod
    def setA2(val):
        LArTrans.A2 = val

    @staticmethod
    def setmu2(val):
        LArTrans.mu2 = val

    @staticmethod
    def setsigma2(val):
        LArTrans.sigma2 = val









