import numpy as np

class LArRindex(object):

    wavelength = [0.3612, 0.3650, 0.4063, 0.4358, 0.4753, 0.5086, 0.5461, 0.5780, 0.6439]  # um
    rindex1    = [1.2395, 1.2392, 1.2372, 1.2361, 1.2349, 1.2341, 1.2334, 1.2328, 1.2321]  # 83.81 K
    rindex2    = [1.2370, 1.2367, 1.2347, 1.2336, 1.2324, 1.2316, 1.2308, 1.2303, 1.2296]  # 86 K 
    rindex3    = [1.2349, 1.2346, 1.2326, 1.2315, 1.2303, 1.2295, 1.2287, 1.2282, 1.2274]  # 88 K
    rindex4    = [1.2326, 1.2331, 1.2308, 1.2297, 1.2285, 1.2277, 1.2269, 1.2264, 1.2256]  # 90 K

    lUV = 0.1066
    lIR = 0.9083
    a0 = 0.335
    aUV = 0.099
    aIR = 0.008
    rho90K = 0.03449
    rhop0 = -0.000158
    rhop1 = 0.0487
    T = 86

    a00 = 0.335
    aUV0 = 0.099
    aIR0 = 0.008
    e_a0 = 0.003
    e_aUV = 0.003
    e_aIR = 0.003
    rhop00 = -0.000158
    rhop10 = 0.0487
    e_rhop0 = 0.000004
    e_rhop1 = 0.0003

    rindex_data = []
    rindex_err  = []
    rindex_calc = []
    
    rindex_data_T1, rindex_err_T1, rindex_calc_T1 = [], [], []
    rindex_data_T2, rindex_err_T2, rindex_calc_T2 = [], [], []
    rindex_data_T3, rindex_err_T3, rindex_calc_T3 = [], [], []
    rindex_data_T4, rindex_err_T4, rindex_calc_T4 = [], [], []
    rho1 = 0.03549    # 83.81 K
    rho2 = 0.03513    # 86 K
    rho3 = 0.03481    # 88 K
    rho4 = 0.03449    # 90 K
    

    chi2 = 0

    @staticmethod
    def LoadData():
        LArRindex.rindex_data = []
        LArRindex.rindex_err  = []
        for i in range(len(LArRindex.wavelength)):
            LArRindex.rindex_data.append((LArRindex.rindex1[i]+LArRindex.rindex2[i]+LArRindex.rindex3[i]) / 3.)
            d1 = abs(LArRindex.rindex1[i]-LArRindex.rindex_data[i])
            d2 = abs(LArRindex.rindex_data[i]-LArRindex.rindex3[i])
            if d1 > d2:
                LArRindex.rindex_err.append(d1)
            else:
                LArRindex.rindex_err.append(d2)

        LArRindex.rindex_data_T1 = []
        LArRindex.rindex_err_T1  = []
        LArRindex.rindex_data_T2 = []
        LArRindex.rindex_err_T2  = []
        LArRindex.rindex_data_T3 = []
        LArRindex.rindex_err_T3  = []
        LArRindex.rindex_data_T4 = []
        LArRindex.rindex_err_T4  = []
        for i in range(len(LArRindex.wavelength)):
            LArRindex.rindex_data_T1.append(LArRindex.rindex1[i])
            LArRindex.rindex_err_T1.append(LArRindex.rindex1[i]*0.001)
            LArRindex.rindex_data_T2.append(LArRindex.rindex2[i])
            LArRindex.rindex_err_T2.append(LArRindex.rindex2[i]*0.001)
            LArRindex.rindex_data_T3.append(LArRindex.rindex3[i])
            LArRindex.rindex_err_T3.append(LArRindex.rindex3[i]*0.001)
            LArRindex.rindex_data_T4.append(LArRindex.rindex4[i])
            LArRindex.rindex_err_T4.append(LArRindex.rindex4[i]*0.001)
    


    @staticmethod
    def rindex_func(l):
        lUV = LArRindex.lUV
        lIR = LArRindex.lIR
        a0 = LArRindex.a0
        aUV = LArRindex.aUV
        aIR = LArRindex.aIR
        T = LArRindex.T
        p0 = LArRindex.rhop0
        p1 = LArRindex.rhop1
        rho = (p0*T+p1) / LArRindex.rho90K
        A = (a0 + aUV*l**2/(l**2-lUV**2) + aIR*l**2/(l**2-lIR**2)) * rho
        n = np.sqrt(1+3*A/(3-A))
        return n 


    @staticmethod
    def rindex_func_density(l, rho):
        lUV = LArRindex.lUV
        lIR = LArRindex.lIR
        a0 = LArRindex.a0
        aUV = LArRindex.aUV
        aIR = LArRindex.aIR
        T = LArRindex.T
        A = (a0 + aUV*l**2/(l**2-lUV**2) + aIR*l**2/(l**2-lIR**2)) * rho
        n = np.sqrt(1+3*A/(3-A))
        return n 



    @staticmethod
    def Calculate():
        LArRindex.rindex_calc = []
        for i in LArRindex.wavelength:
            LArRindex.rindex_calc.append(LArRindex.rindex_func(i))

    @staticmethod
    def Calculate_new():
        LArRindex.rindex_calc_T1 = []
        LArRindex.rindex_calc_T2 = []
        LArRindex.rindex_calc_T3 = []
        LArRindex.rindex_calc_T4 = []
        
        for i in LArRindex.wavelength:
            LArRindex.rindex_calc_T1.append(LArRindex.rindex_func_density(i, LArRindex.rho1/LArRindex.rho90K))
            LArRindex.rindex_calc_T2.append(LArRindex.rindex_func_density(i, LArRindex.rho2/LArRindex.rho90K))
            LArRindex.rindex_calc_T3.append(LArRindex.rindex_func_density(i, LArRindex.rho3/LArRindex.rho90K))
            LArRindex.rindex_calc_T4.append(LArRindex.rindex_func_density(i, LArRindex.rho4/LArRindex.rho90K))


    
    @staticmethod
    def GetChi2():
        LArRindex.chi2 = 0
        LArRindex.Calculate()
        for i in range(len(LArRindex.rindex_data)):
            LArRindex.chi2 += (LArRindex.rindex_data[i]-LArRindex.rindex_calc[i])**2/LArRindex.rindex_err[i]**2
    
        return LArRindex.chi2


    @staticmethod 
    def GetChi2_new():
        LArRindex.chi2 = 0
        LArRindex.Calculate_new()
        for i in range(len(LArRindex.wavelength)):
            LArRindex.chi2 += (LArRindex.rindex_data_T1[i]-LArRindex.rindex_calc_T1[i])**2 / LArRindex.rindex_err_T1[i]**2
            LArRindex.chi2 += (LArRindex.rindex_data_T2[i]-LArRindex.rindex_calc_T2[i])**2 / LArRindex.rindex_err_T2[i]**2
            LArRindex.chi2 += (LArRindex.rindex_data_T3[i]-LArRindex.rindex_calc_T3[i])**2 / LArRindex.rindex_err_T3[i]**2
            LArRindex.chi2 += (LArRindex.rindex_data_T4[i]-LArRindex.rindex_calc_T4[i])**2 / LArRindex.rindex_err_T4[i]**2
        return LArRindex.chi2
    

    @staticmethod
    def GetPulls():
        pull = 0
        pull += (LArRindex.a0 - LArRindex.a00)**2 / LArRindex.e_a0**2
        pull += (LArRindex.aUV - LArRindex.aUV0)**2 / LArRindex.e_aUV**2
        pull += (LArRindex.aIR - LArRindex.aIR0)**2 / LArRindex.e_aIR**2
        pull += (LArRindex.rhop0 - LArRindex.rhop00)**2 / LArRindex.e_rhop0**2
        pull += (LArRindex.rhop1 - LArRindex.rhop10)**2 / LArRindex.e_rhop1**2
        return pull


    @staticmethod
    def Plot():
        print("-------------------> Saving fitted refractive index values ! +++++++++++++ ")
        import matplotlib.pyplot as plt
        from matplotlib import gridspec
        LArRindex.Calculate_new()
        fig = plt.figure(figsize=(11 ,9))
        spec = gridspec.GridSpec(ncols=2, nrows=2)
        ax1 = fig.add_subplot(spec[0])
        ax2 = fig.add_subplot(spec[1])
        ax3 = fig.add_subplot(spec[2])
        ax4 = fig.add_subplot(spec[3])
        ax1.errorbar(LArRindex.wavelength, LArRindex.rindex_data_T1, yerr=LArRindex.rindex_err_T1, fmt="o", color="blue", label="Simulation: T=83.81 K")
        ax1.plot(LArRindex.wavelength, LArRindex.rindex_calc_T1, "v", color="blue", fillstyle="none", label="Calculation: T=83.81 K")
        ax2.errorbar(LArRindex.wavelength, LArRindex.rindex_data_T2, yerr=LArRindex.rindex_err_T2, fmt="o", color="coral", label="Simulation: T=86 K")
        ax2.plot(LArRindex.wavelength, LArRindex.rindex_calc_T2, "v", color="coral", fillstyle="none", label="Calculation: T=86 K")
        ax3.errorbar(LArRindex.wavelength, LArRindex.rindex_data_T1, yerr=LArRindex.rindex_err_T3, fmt="o", color="forestgreen", label="Simulation: T=88 K")
        ax3.plot(LArRindex.wavelength, LArRindex.rindex_calc_T3, "v", color="forestgreen", fillstyle="none", label="Calculation: T=88 1K")
        ax4.errorbar(LArRindex.wavelength, LArRindex.rindex_data_T4, yerr=LArRindex.rindex_err_T4, fmt="o", color="violet", label="Simulation: T=90 K")
        ax4.plot(LArRindex.wavelength, LArRindex.rindex_calc_T4, "v", color="violet", fillstyle="none", label="Calculation: T=90 K")
        dx = np.arange(0.120, 0.7, 0.001)
        dy1 = LArRindex.rindex_func_density(dx, LArRindex.rho1/LArRindex.rho90K)
        dy2 = LArRindex.rindex_func_density(dx, LArRindex.rho2/LArRindex.rho90K)
        dy3 = LArRindex.rindex_func_density(dx, LArRindex.rho3/LArRindex.rho90K)
        dy4 = LArRindex.rindex_func_density(dx, LArRindex.rho4/LArRindex.rho90K)
        ax1.plot(dx, dy1, "-", lw=1.5, color="gray")
        ax2.plot(dx, dy2, "-", lw=1.5, color="gray")
        ax3.plot(dx, dy3, "-", lw=1.5, color="gray")
        ax4.plot(dx, dy4, "-", lw=1.5, color="gray")
        ax1.plot(0.128, LArRindex.rindex_func_density(0.128, LArRindex.rho1/LArRindex.rho90K), "*", ms=8, color="red")
        ax2.plot(0.128, LArRindex.rindex_func_density(0.128, LArRindex.rho2/LArRindex.rho90K), "*", ms=8, color="red")
        ax3.plot(0.128, LArRindex.rindex_func_density(0.128, LArRindex.rho3/LArRindex.rho90K), "*", ms=8, color="red")
        ax4.plot(0.128, LArRindex.rindex_func_density(0.128, LArRindex.rho4/LArRindex.rho90K), "*", ms=8, color="red")
        #ax.text(0.131, LArRindex.rindex_func(0.128), "n(128nm)=%.3f"%LArRindex.rindex_func(0.128), fontsize=14)
        ax1.legend(prop={"size":14})
        ax1.set_xlabel("wavelength [um]", fontsize=14)
        ax1.set_ylabel("refractive index", fontsize=14)
        ax1.grid(True)
        ax2.legend(prop={"size":14})
        ax2.set_xlabel("wavelength [um]", fontsize=14)
        ax2.set_ylabel("refractive index", fontsize=14)
        ax2.grid(True)
        ax3.legend(prop={"size":14})
        ax3.set_xlabel("wavelength [um]", fontsize=14)
        ax3.set_ylabel("refractive index", fontsize=14)
        ax3.grid(True)
        ax4.legend(prop={"size":14})
        ax4.set_xlabel("wavelength [um]", fontsize=14)
        ax4.set_ylabel("refractive index", fontsize=14)
        ax4.grid(True)

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
    def getrhop0():
        return LArRindex.rhop0

    @staticmethod
    def getrhop1():
        return LArRindex.rhop1

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

    @staticmethod
    def getrho90K():
        return LArRindex.rho90K

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
    def setrhop0(val):
        LArRindex.rhop0 = val

    @staticmethod
    def setrhop1(val):
        LArRindex.rhop1 = val

    @staticmethod
    def setT(val):
        LArRindex.T = val



