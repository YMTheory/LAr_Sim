#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a*(x-b)**2+c

def find_interval(xarr, yarr):
    for i in range(len(yarr)-1):
        if yarr[i] >= 1 and yarr[i+1] <= 1:
            left = xarr[i] + (yarr[i]-1) / (yarr[i]-yarr[i+1]) * (xarr[i+1]-xarr[i])
            print("Left margin = %.3f" %left)
        if yarr[i] <= 1 and yarr[i+1] >= 1:
            right = xarr[i] + (yarr[i]-1) / (yarr[i]-yarr[i+1]) * (xarr[i+1]-xarr[i])
            print("Right margin = %.3f" %right)


def main():

    chimin = 206.785

    rindex, chi2 = [], []
    with open('../rindex_lag.txt') as f:
        for lines in f.readlines():
            line = lines.strip('\n')
            data = line.split(' ')
            rindex.append(float(data[0]))
            chi2.append(float(data[1])-chimin)

    find_interval(rindex, chi2)

    popt1, pcov1 = curve_fit(func, rindex, chi2, p0=[344860, 60, 208])
    print(popt1)

    #drawx = np.arange(55, 70, 0.001)
    drawx = np.arange(1.340, 1.372, 0.0001)
    drawy = [func(i, popt1[0], popt1[1], popt1[2]) for i in drawx]

    plt.plot(drawx, drawy, "--", lw=2)

    plt.plot(rindex, chi2, "o-", color='royalblue', ms=2)
    plt.text(1.3555, 210.0, r"$344864.229\times(x-1.357)^2+208.750$", fontsize=13)
    #plt.text(56, 210.2, r"$0.0389\times(x-61.568)^2+208.70$", fontsize=13)
    plt.hlines(1, 50, 80, color='darkviolet')
    #plt.hlines(popt1[2]+1, 1.355, 1.372, color='darkviolet')

    plt.xlabel(r"$L_{Ray} /m$")
    #plt.xlabel("rindex")
    plt.ylabel(r"$\chi^2$")

    plt.xlim(50, 80)
    #plt.xlim(1.357, 1.372)
    plt.ylim(-1, 6)

    plt.tight_layout()
    plt.grid(True)

    #plt.savefig("Lagrange_rayL128_new.pdf")
    plt.show()


if __name__ == '__main__':
    plt.style.use("seaborn-muted")
    main()
