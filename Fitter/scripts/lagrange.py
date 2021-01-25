#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a*(x-b)**2+c

def main():

    chimin = 208.343

    rindex, chi2 = [], []
    with open('./rindex_1dscan.txt') as f:
        for lines in f.readlines():
            line = lines.strip('\n')
            data = line.split(' ')
            rindex.append(float(data[0]))
            chi2.append(float(data[1])-chimin)

    popt1, pcov1 = curve_fit(func, rindex, chi2, p0=[344860, 1.357, 208])
    print(popt1)

    #drawx = np.arange(55, 70, 0.001)
    drawx = np.arange(1.340, 1.370, 0.0001)
    drawy = [func(i, popt1[0], popt1[1], popt1[2]) for i in drawx]

    plt.plot(drawx, drawy, "--", lw=2)

    plt.plot(rindex, chi2, "o", ms=2)
    plt.text(1.3555, 210.0, r"$344864.229\times(x-1.357)^2+208.750$", fontsize=13)
    #plt.text(56, 210.2, r"$0.0389\times(x-61.568)^2+208.70$", fontsize=13)
    plt.hlines(popt1[2]+1, 1.355, 1.359, color='darkviolet')

    plt.xlabel("rindex")
    plt.ylabel(r"$\chi^2$")

    plt.xlim(1.346, 1.368)
    plt.ylim(-1, 10)

    #plt.savefig("Lagrange_rindex128.png")
    plt.show()


if __name__ == '__main__':
    plt.style.use("seaborn-muted")
    main()