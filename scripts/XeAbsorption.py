#!/usr/bin/env python
# coding=utf-8

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def readData(filename):
    wavelength, absorbance = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            wavelength.append(float(data[0]))
            absorbance.append(float(data[1]))
    return wavelength, absorbance

def gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)*(x-mu)/2/sigma/sigma)

if __name__ == '__main__':
    wavelength, absorbance = readData("XeDopedAbs.csv")

    # locallly peak fitting :
    wl1peak, abs1peak, wl2peak, abs2peak = [], [], [], []
    peak1, peak2, width1, width2 = 127.0, 141.5, 1.5, 5.
    for wl, ab in zip(wavelength, absorbance):
        if abs(wl-peak1)<=width1:
            wl1peak.append(wl)
            abs1peak.append(ab)
        elif abs(wl-peak2)<=width2:
            wl2peak.append(wl)
            abs2peak.append(ab)

    popt1, pcov1 = curve_fit(gauss, wl1peak, abs1peak, p0=[1, peak1, 1])
    popt2, pcov2 = curve_fit(gauss, wl2peak, abs2peak, p0=[1, peak2, 1])

    print(popt2)
    print(pcov2)

    drawx1, drawy1, drawx2, drawy2 = [], [], [], []
    for i in range(100):
        drawx1.append(peak1-width1+width1*2/100.*i)
        drawy1.append(gauss(drawx1[-1], popt1[0], popt1[1], popt1[2]))
        drawx2.append(peak2-width2+width2*2/100.*i)
        drawy2.append(gauss(drawx2[-1], popt2[0], popt2[1], popt2[2]))

    plt.plot(wavelength, absorbance, "o--", ms=5, color="blue")
    plt.plot(drawx1, drawy1, "-", color="green")
    plt.text(120, 1.1, r"$ %.3f e^{-\frac{(x-%.3f)^2}{2 \times %.2f ^2}}$" %(popt1[0], popt1[1], popt1[2]), fontsize=14, color='green')
    plt.plot(drawx2, drawy2, "-", color="orange")
    plt.text(130, 1.0, r"$ %.3f e^{-\frac{(x-%.3f)^2}{2 \times %.2f ^2}}$" %(popt2[0], popt2[1], popt2[2]), fontsize=14, color='orange')
    plt.xlabel("wavelenght/nm")
    plt.ylabel("absorbance")
    #plt.savefig("XeDopedAbs.png")
    #plt.show()