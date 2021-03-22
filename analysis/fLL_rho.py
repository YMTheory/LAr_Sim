#!/usr/bin/env python
# coding=utf-8

# liquid 
T1, rho1 = 83.81, 3.549
T2, rho2 = 86, 3.513
T3, rho3 = 88, 3.481
T4, rho4 = 90, 3.449

T = [T1, T2, T3, T4]
rho = [rho1, rho2, rho3, rho4]

wl = [361.2, 365, 406.3, 435.8, 475.3, 508.6, 546.1, 578, 643.9]
rindex1 = [1.2395, 1.2392, 1.2372, 1.2361, 1.2349, 1.2341, 1.2334, 1.2328, 1.2321]
rindex2 = [1.2370, 1.2367, 1.2347, 1.2336, 1.2324, 1.2316, 1.2308, 1.2303, 1.2295]
rindex3 = [1.2349, 1.2346, 1.2326, 1.2315, 1.2303, 1.2295, 1.2287, 1.2282, 1.2274]
rindex4 = [1.2326, 1.2331, 1.2308, 1.2297, 1.2285, 1.2277, 1.2269, 1.2264, 1.2256]

import numpy as np
rindex1 = np.array(rindex1)
rindex2 = np.array(rindex2)
rindex3 = np.array(rindex3)
rindex4 = np.array(rindex4)


def draw(x, y, T, rho):
    plt.plot(x, y, "o-", label=r"$T=%.2f K, \rho=%.3f 10^{-2}mole/cm^3 $" %(T, rho))

def calc_fLL():
    fLL1 = (rindex1**2 - 1)/(rindex1**2 +2) * 1/rho1
    fLL2 = (rindex2**2 - 1)/(rindex2**2 +2) * 1/rho2
    fLL3 = (rindex3**2 - 1)/(rindex3**2 +2) * 1/rho3
    fLL4 = (rindex4**2 - 1)/(rindex4**2 +2) * 1/rho4
    return fLL1, fLL2, fLL3, fLL4

import matplotlib.pyplot as plt
def draw_rindex():
    plt.style.use("seaborn-deep")
    draw(wl, rindex1, T1, rho1)
    draw(wl, rindex2, T2, rho2)
    draw(wl, rindex3, T3, rho3)
    draw(wl, rindex4, T4, rho4)

    plt.ylabel("refractive index")
    plt.xlabel("wavelength/nm")
    plt.title("Sinnock 1969")
    plt.legend()
    plt.grid(True)
    #plt.savefig("rindex_compare.pdf")
    plt.show()

def draw_fLL():
    plt.style.use("seaborn-deep")
    fLL1, fLL2, fLL3, fLL4 = calc_fLL()

    fLL_wl0, fLL_wl1, fLL_wl2, fLL_wl3 = [], [], [], []

    fLL_wl0.append(fLL1[0])
    fLL_wl0.append(fLL2[0])
    fLL_wl0.append(fLL3[0])
    fLL_wl0.append(fLL4[0])
    fLL_wl1.append(fLL1[2])
    fLL_wl1.append(fLL2[2])
    fLL_wl1.append(fLL3[2])
    fLL_wl1.append(fLL4[2])
    fLL_wl2.append(fLL1[4])
    fLL_wl2.append(fLL2[4])
    fLL_wl2.append(fLL3[4])
    fLL_wl2.append(fLL4[4])
    fLL_wl3.append(fLL1[6])
    fLL_wl3.append(fLL2[6])
    fLL_wl3.append(fLL3[6])
    fLL_wl3.append(fLL4[6])

    mean0 = np.array(fLL_wl0).mean()
    mean1 = np.array(fLL_wl1).mean()
    mean2 = np.array(fLL_wl2).mean()
    mean3 = np.array(fLL_wl3).mean()

    #draw(wl, fLL1, T1, rho1)
    #draw(wl, fLL2, T2, rho2)
    #draw(wl, fLL3, T3, rho3)
    #draw(wl, fLL4, T4, rho4)
    #plt.xlabel("wavelength/nm")
    #plt.title("Sinnock 1969")

    plt.plot(rho, fLL_wl0, "o-", label="%.1f nm"%wl[0])
    plt.fill_between(rho, [mean0*0.999 for i in range(4)], [mean0*1.001 for i in range(4)], alpha=0.3)
    plt.plot(rho, fLL_wl1, "o-", label="%.1f nm"%wl[2])
    plt.fill_between(rho, [mean1*0.998 for i in range(4)], [mean1*1.002 for i in range(4)], alpha=0.3)
    plt.plot(rho, fLL_wl2, "o-", label="%.1f nm"%wl[4])
    plt.fill_between(rho, [mean2*0.998 for i in range(4)], [mean2*1.002 for i in range(4)], alpha=0.3)
    plt.plot(rho, fLL_wl3, "o-", label="%.1f nm"%wl[6])
    plt.fill_between(rho, [mean3*0.998 for i in range(4)], [mean3*1.002 for i in range(4)], alpha=0.3)
    #plt.xlabel("temperature/K")
    plt.xlabel(r"density($10^{-2} mole/cm^3$)")

    plt.ylabel(r"$f_{LL}$")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("fLL2density_compare.pdf")
    #plt.show()

def rindex2_rho():
    plt.style.use("seaborn-deep")

    rho_array = [rho1, rho2, rho3, rho4]
    r0 = [rindex1[0]**2, rindex2[0]**2, rindex3[0]**2, rindex4[0]**2]
    r1 = [rindex1[1]**2, rindex2[1]**2, rindex3[1]**2, rindex4[1]**2]
    r2 = [rindex1[2]**2, rindex2[2]**2, rindex3[2]**2, rindex4[2]**2]
    r3 = [rindex1[3]**2, rindex2[3]**2, rindex3[3]**2, rindex4[3]**2]
    r4 = [rindex1[4]**2, rindex2[4]**2, rindex3[4]**2, rindex4[4]**2]
    r5 = [rindex1[5]**2, rindex2[5]**2, rindex3[5]**2, rindex4[5]**2]
    r6 = [rindex1[6]**2, rindex2[6]**2, rindex3[6]**2, rindex4[6]**2]
    r7 = [rindex1[7]**2, rindex2[7]**2, rindex3[7]**2, rindex4[7]**2]
    r8 = [rindex1[8]**2, rindex2[8]**2, rindex3[8]**2, rindex4[8]**2]

    plt.plot(rho_array, r0, "o-", label="wl=%.1f nm"%wl[0])
    #plt.plot(rho_array, r1, "o-", label="wl=%.1f nm"%wl[1])
    #plt.plot(rho_array, r2, "o-", label="wl=%.1f nm"%wl[2])
    #plt.plot(rho_array, r3, "o-", label="wl=%.1f nm"%wl[3])
    #plt.plot(rho_array, r4, "o-", label="wl=%.1f nm"%wl[4])
    #plt.plot(rho_array, r5, "o-", label="wl=%.1f nm"%wl[5])
    #plt.plot(rho_array, r6, "o-", label="wl=%.1f nm"%wl[6])
    #plt.plot(rho_array, r7, "o-", label="wl=%.1f nm"%wl[7])
    #plt.plot(rho_array, r8, "o-", label="wl=%.1f nm"%wl[8])

    from scipy.optimize import curve_fit
    def func(x, a, b):
        return a*x+b

    popt, pcov = curve_fit(func, rho_array, r0)
    print(popt)
    drawx, drawy = np.arange(rho4, rho1, 0.001), []
    for i in drawx:
        drawy.append(func(i, *popt))
    plt.plot(drawx, drawy, "--")

    plt.xlabel(r"$\rho 10^{-2} mole/cm^3$")
    plt.ylabel(r"n^2")
    plt.grid(True)
    plt.legend()
    plt.show()


def main():
    draw_fLL()
    #rindex2_rho()

if __name__ == '__main__':
    main()