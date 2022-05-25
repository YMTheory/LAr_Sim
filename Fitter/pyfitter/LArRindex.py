import numpy as np

class LArRindex(object):

    wavelength = [0.3612, 0.3650, 0.4063, 0.4358, 0.4753, 0.5086, 0.5461, 0.5780, 0.6439]  # um
    rindex0    = [1.2395, 1.2392, 1.2372, 1.2361, 1.2349, 1.2341, 1.2334, 1.2328, 1.2321]  # 83.81 K
    rindex1    = [1.2370, 1.2367, 1.2347, 1.2336, 1.2324, 1.2316, 1.2308, 1.2303, 1.2296]  # 86 K 
    rindex2    = [1.2349, 1.2346, 1.2326, 1.2315, 1.2303, 1.2295, 1.2287, 1.2282, 1.2274]  # 88 K

    lUV = 0.1066
    lIR = 0.9083
    a0 = 0.335
    aUV = 0.099
    aIR = 0.008
    rho90K = 0.03449
    rho0 = -0.000158
    rho1 = 0.0487
    T = 86

    rindex_data = []
    rindex_err  = []
    rindex_calc = []
    chi2 = 0

    @staticmethod
    def LoadData():
        LArRindex.rindex_data = []
        LArRindex.rindex_err  = []
        for i in range(len(LArRindex.wavelength)):
            LArRindex.rindex_data.append((LArRindex.rindex0[i]+LArRindex.rindex1[i]+LArRindex.rindex2[i]) / 3.)
            d1 = abs(LArRindex.rindex0[i]-LArRindex.rindex_data[i])
            d2 = abs(LArRindex.rindex_data[i]-LArRindex.rindex2[i])
            if d1 > d2:
                LArRindex.rindex_err.append(d1)
            else:
                LArRindex.rindex_err.append(d2)



    @staticmethod
    def rindex_func(l):
        lUV = LArRindex.lUV
        lIR = LArRindex.lIR
        a0 = LArRindex.a0
        aUV = LArRindex.aUV
        aIR = LArRindex.aIR
        T = LArRindex.T
        p0 = LArRindex.rho0
        p1 = LArRindex.rho1
        rho = (p0*T+p1) / LArRindex.rho90K
        A = (a0 + aUV*l**2/(l**2-lUV**2) + aIR*l**2/(l**2-lIR**2)) * rho
        n = np.sqrt(1+3*A/(3-A))
        return n 


    @staticmethod
    def Calculate():
        LArRindex.rindex_calc = []
        for i in LArRindex.wavelength:
            LArRindex.rindex_calc.append(LArRindex.rindex_func(i))
        
    
    @staticmethod
    def GetChi2():
        LArRindex.chi2 = 0
        LArRindex.Calculate()
        for i in range(len(LArRindex.rindex_data)):
            LArRindex.chi2 += (LArRindex.rindex_data[i]-LArRindex.rindex_calc[i])**2/LArRindex.rindex_err[i]**2
    
        return LArRindex.chi2


    @staticmethod
    def Plot():
        import matplotlib.pyplot as plt
        LArRindex.Calculate()
        fig, ax = plt.subplots()
        ax.errorbar(LArRindex.wavelength, LArRindex.rindex_data, yerr=LArRindex.rindex_err, fmt="o", label="Simulation")
        ax.plot(LArRindex.wavelength, LArRindex.rindex_calc, "o", label="Calculation")
        dx = np.arange(0.120, 0.7, 0.001)
        dy = LArRindex.rindex_func(dx)
        ax.plot(dx, dy, "-", lw=1.5, color="gray")
        ax.plot(0.128, LArRindex.rindex_func(0.128), "*", ms=8, color="red")
        ax.text(0.131, LArRindex.rindex_func(0.128), "n(128nm)=%.3f"%LArRindex.rindex_func(0.128), fontsize=14)
        ax.legend(prop={"size":14})
        ax.set_xlabel("wavelength [um]", fontsize=14)
        ax.set_ylabel("refractive index", fontsize=14)
        ax.grid(True)
        plt.tight_layout()
        plt.savefig("rindex.pdf")
                

    ######## Getter Functions ######
    @staticmethod
    def getlUV():
        return LArRindex.lUV
    
    @staticmethod
    def getlIR():
        return LArRindex.lIR

    @staticmethod
    def geta0():
        return LArRindex.a0

    @staticmethod
    def getaUV():
        return LArRindex.aUV
    
    @staticmethod
    def getaIR():
        return LArRindex.aIR

    @staticmethod
    def getrho0():
        return LArRindex.rho0

    @staticmethod
    def getrho1():
        return LArRindex.rho1

    @staticmethod
    def getT():
        return LArRindex.T

    @staticmethod
    def getDataX():
        return LArRindex.wavelength

    @staticmethod
    def getDataY():
        return LArRindex.rindex_data

    @staticmethod
    def getDataYerr():
        return LArRindex.rindex_err

    ##########  Setter Functions ##########
    @staticmethod
    def setlUV(val):
        LArRindex.lUV = val

    @staticmethod
    def setlIR(val):
        LArRindex.lIR = val

    @staticmethod
    def seta0(val):
        LArRindex.a0 = val

    @staticmethod
    def setaUV(val):
        LArRindex.aUV = val

    @staticmethod
    def setaIR(val):
        LArRindex.aIR = val

    @staticmethod
    def setrho0(val):
        LArRindex.rho0 = val

    @staticmethod
    def setrho1(val):
        LArRindex.rho1 = val

    @staticmethod
    def setT(val):
        LArRindex.T = val



