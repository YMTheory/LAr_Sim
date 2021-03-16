#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

def load(filename):
    pressure, volume = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            pressure.append(float(data[1]))
            volume.append(float(data[3]))
    
    return pressure, volume

T = [84, 85, 86 ,87, 90] # 5 conditions in 1 atm

def linearfunc(x, a, b):
    return a*x + b

from scipy.optimize import curve_fit

def main():
    plt.style.use("seaborn-muted")

    x = np.arange(1.32, 1.5, 0.01)

    """
    p0, v0 = load("data0.txt")
    popt0, pcov0 = curve_fit(linearfunc, p0, v0)
    y0 = linearfunc(x, *popt0)
    plt.plot(p0, v0, "o", label="T=%d K"%T[0])
    plt.plot(x, y0, "-", color="royalblue")
    """
    """
    p1, v1 = load("data1.txt")
    popt1, pcov1 = curve_fit(linearfunc, p1, v1)
    y1 = linearfunc(x, *popt1)
    plt.plot(p1, v1, "o", label="T=%d K"%T[1], color="limegreen")
    plt.plot(x, y1, "-", color="limegreen")
    """

    """
    p2, v2 = load("data2.txt")
    popt2, pcov2 = curve_fit(linearfunc, p2, v2)
    y2 = linearfunc(x, *popt2)
    plt.plot(p2, v2, "o", label="T=%d K"%T[2], color="peru")
    plt.plot(x, y2, color="peru")

    p3, v3 = load("data3.txt")
    popt3, pcov3 = curve_fit(linearfunc, p3, v3)
    y3 = linearfunc(x, *popt3)
    plt.plot(p3, v3, "o", label="T=%d K"%T[3], color="darkviolet")
    plt.plot(x, y3, color="darkviolet")
    """
    
    """

    p4, v4 = load("data4.txt")
    popt4, pcov4 = curve_fit(linearfunc, p4, v4)
    y4 = linearfunc(x, *popt4)
    plt.plot(p4, v4, "o", label="T=%d K"%T[4], color="darkviolet")
    plt.plot(x, y4, color="darkviolet")

    plt.text(x[4], y4[2], "slope=%e" %popt4[0], fontsize=12)
    """

    kappaT = [1.9355e-9, 1.9903e-9, 2.0475e-9, 2.1064e-9, 2.2983e-9]  # m2/N

    popt, pcov = curve_fit(linearfunc, T, kappaT)
    print(popt)
    print(pcov)
    dx = np.arange(83, 91, 0.1)
    dy = linearfunc(dx, *popt)
    plt.plot(dx, dy, "--")

    plt.plot(T, kappaT, "o")
    plt.xlabel("temperature/K")
    plt.ylabel(r"$m^2/N$")    
    plt.title("LAr isothermal compressibility")

    plt.text(T[0], kappaT[4], r"$\kappa_T$ = %e*T+%e" %(popt[0], popt[1]), fontsize=14)

    #plt.xlabel("pressure/atm")
    #plt.ylabel("L/mol")

    plt.tight_layout()
    #plt.legend()
    plt.grid(True)
    #ax.set_axisbelow(True)

    plt.savefig("kappaT_T.pdf")
    plt.show()


if __name__=="__main__":
    main()
