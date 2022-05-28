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

    R = 0.937
    A1  = 0.4
    mu1 = 126.51
    sigma1 = 1
    A2 = 0.4
    mu2 = 140.12
    sigma2 = 1.537

    nu_f = 0.0
    e_f = 1
    k00 = 6.0712e-11
    k10 = -3.1699e-09
    e_k0 = np.sqrt(1.49189964e-24)
    e_k1 = np.sqrt(1.11432969e-20)
    T0 = 86
    e_T = 3./np.sqrt(12.)

    chi2 = 0

    @staticmethod
    def LoadData():
        LArTrans.wavelength = []
        LArTrans.trans_data = []
        LArTrans.trans_err  = []
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
        LArRindex.setT(LArTrans.T)
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
        LArRindex.setT(LArTrans.T)
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
        #A1     = LArTrans.A1
        R      = LArTrans.R
        mu1    = LArTrans.mu1
        sigma1 = LArTrans.sigma1
        A2     = LArTrans.A2
        mu2    = LArTrans.mu2
        sigma2 = LArTrans.sigma2
        A1 = R * A2

        A_abs = A1*np.exp(-(l-mu1)*(l-mu1)/2/sigma1/sigma1) + A2*np.exp(-(l-mu2)*(l-mu2)/2/sigma2/sigma2);
        T_abs = np.exp( -A_abs*np.log(10.) );
        return T_abs

    @staticmethod
    def fresnel_func(l):
        wl = l * 1000
        E = 1240./wl
        # MgF2, -193 degree
        a = 29
        E0 = 12.0226
        gamma = 0.3594
        b = 295.855
        E1 = 20.8448
        
        n_MgF2 = np.sqrt(1+a*(E0*E0-E*E)/((E0*E0-E*E)*(E0*E0-E*E)+gamma*E*E) + b/(E1*E1-E*E))

        n_vac = 1.
        LArRindex.setT(LArTrans.T)
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
        LArRindex.setT(LArTrans.T)
        LArTrans.trans_calc = []
        for i in LArTrans.wavelength:
            trans_pred = LArTrans.transmission_func(i)
            LArTrans.trans_calc.append(trans_pred * (1 + LArTrans.nu_f))


    @staticmethod
    def GetChi2():
        LArTrans.chi2 = 0
        LArTrans.Calculate()
        for i in range(len(LArTrans.wavelength)):
            LArTrans.chi2 += (LArTrans.trans_data[i] - LArTrans.trans_calc[i])**2 / (LArTrans.trans_err[i])**2
        return LArTrans.chi2


    @staticmethod
    def GetPulls():
        pull = 0
        pull += (LArTrans.k0 - LArTrans.k00)**2 / LArTrans.e_k0**2
        pull += (LArTrans.k1 - LArTrans.k10)**2 / LArTrans.e_k1**2
        pull += LArTrans.nu_f**2 / LArTrans.e_f**2
        pull += (LArTrans.T - LArTrans.T0)**2 / LArTrans.e_T**2

        return pull



    @staticmethod
    def Plot():
        print("-------------------> Saving fitted transmission spectrum ! +++++++++++++ ")
        import matplotlib.pyplot as plt
        LArTrans.Calculate()
        fig, ax = plt.subplots()
        draw_data, draw_err = [], []
        for i in range(len(LArTrans.wavelength)):
            draw_data.append(LArTrans.trans_data[i] / LArTrans.fresnel_func(LArTrans.wavelength[i]))
            draw_err.append(LArTrans.trans_err[i] / LArTrans.fresnel_func(LArTrans.wavelength[i]))
        ax.errorbar(LArTrans.wavelength, draw_data, yerr=draw_err, fmt="o", label="Simulation")
        dx = np.arange(0.125, 0.150, 0.0001)
        draw_calc = []
        for i in range(len(LArTrans.wavelength)):
            draw_calc.append(LArTrans.Tabs_func(LArTrans.wavelength[i]) * np.exp(-5.8/LArTrans.lray_func(LArTrans.wavelength[i])))
        ax.plot(LArTrans.wavelength, draw_calc, "^", fillstyle="none", label="Calc: tot")
        ax.plot(dx, np.exp(-5.8/LArTrans.lray_func(dx)), "-", label="Calc: Ray")
        ax.plot(dx, LArTrans.Tabs_func(dx), "-", label="Calc: Abs")
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

    @staticmethod
    def getnuf():
        return LArTrans.nu_f

    @staticmethod
    def getR():
        return LArTrans.R


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

    @staticmethod
    def setnuf(val):
        LArTrans.nu_f = val


    @staticmethod
    def setR(val):
        LArTrans.R = val





