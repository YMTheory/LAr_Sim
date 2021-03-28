import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn-muted')

def babizcmodel(a0, aUV, aIR, ratio, x):
    A = a0 + aUV*x**2/(x**2-0.1066**2) + aIR*x**2/(x**2-0.9083**2)
    A *= ratio
    return np.sqrt(1+3*A/(3-A))


def finiteDiv(n, l):
    x = 17.5
    y = 11.6
    z = 36.5
    return ((x+y+z)/(x+y/n+z))**2


def readCalc():
    wl_arr, corr_arr = [], []
    with open("./calcFiniteDiv.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            wl_arr.append(float(data[0]))
            corr_arr.append(float(data[1]))
    
    return wl_arr, corr_arr


def FresnelCorr(l, rindex):
    l = l*1000. #nm
    n_vac = 1; # vacuum
    a = 38 #eV2
    E0 = 12.8234 #eV
    gamma = 0.42357 #eV
    b = 247.188 #eV2
    E1 = 18.8448 #eV

    n_MgF2 = np.sqrt(1+a*(E0**2-(1240./l)**2)/((E0**2-(1240/l)**2)**2 + gamma**2*(1240./l)**2) + b/(E1**2-(1240/l)**2))

    factor = ( (1-((n_vac-n_MgF2)/(n_vac+n_MgF2))**2)/(1-((n_MgF2-rindex)/(n_MgF2+rindex))**2) )**2

    return factor


def main():
    corr_arr, wl_arr = [], []
    for wl in range(120, 150, 1):
        wl = wl/1000 #unit um
        wl_arr.append(wl)
        rindex = babizcmodel(0.3347, 0.0994, 0.0080, 0.03507/0.03449, wl)
        corr_arr.append( FresnelCorr(wl, rindex) * finiteDiv(rindex, wl) )

    wl_arr, corr_arr = readCalc()

    plt.plot(wl_arr, corr_arr, "-", color='red')
    
    plt.xlabel("wavelength/nm")
    plt.ylabel("correction factor")

    plt.grid(True)
    #plt.savefig("finiteDiv.pdf")
    plt.show()

if __name__ == '__main__':
    main()