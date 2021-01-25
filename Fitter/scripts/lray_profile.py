#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import math

def load(filename, chimin):
    ratio, lray, chi2 = [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip('\n')
            data = line.split(' ')
            ratio.append(float(data[0]))
            lray.append(float(data[1]))
            chi2.append(float(data[2])-chimin)
    return np.array(ratio), np.array(lray), np.array(chi2)

def cross(a, b, flag):
    if flag == 0:
        return b-np.sqrt(1/a)
    else:
        return b+np.sqrt(1/a)

if __name__ == '__main__':
    plt.style.use("seaborn-deep")

    ratio, lray, chi2 = load('../log2d.txt', 208.397)
    x, y =[], []
    for i in range(len(ratio)):
        if abs(chi2[i] - 4 ) <= 0.03:
            x.append(lray[i])
            y.append(ratio[i])
    print("final contour point number: %d" %len(x) )
    plt.plot(lray, ratio, "s", ms=2)
    plt.plot(x, y, "o", ms=4)
    plt.xlim(40, 80)
    plt.show()


    """
    ratio1, lray1, chi1 = load("../ba.log", 208.397)
    ratio2, lray2, chi2 = load("../our.log", 209.417)

    #x, y = [], []
    #for i in seed:
    #    x.append(lray[int(i)])
    #    y.append(chi[int(i)]-208.397)

    from scipy.optimize import curve_fit
    def func(x, a, b, c):
        return a*(x-b)**2+c

    popt1, pcov1 = curve_fit(func, lray1, chi1, p0=[0.05, 60, 0])    
    print(popt1)
    cross11 = cross(popt1[0], popt1[1], 0)
    cross12 = cross(popt1[0], popt1[1], 1)

    dx1, dy1 = [], []
    for i in range(1000):
        dx1.append(50 + 20./1000*i)
        dy1.append(func(dx1[-1], popt1[0], popt1[1], 0))

    popt2, pcov2 = curve_fit(func, lray2, chi2, p0=[0.08, 62, 0])    
    print(popt2)
    cross21 = cross(popt2[0], popt2[1], 0)
    cross22 = cross(popt2[0], popt2[1], 1)

    dx2, dy2 = [], []
    for i in range(1000):
        dx2.append(50 + 20./1000*i)
        dy2.append(func(dx2[-1], popt2[0], popt2[1], 0))

    #plt.plot(lray1, chi1-208.397, "o", ms=4)
    plt.plot([60.139, 60.139], [-2, 10], "--", color='blue', label="model1 best-fit")
    plt.plot(dx1, dy1, '-')
    plt.fill_betweenx([0, 8], [cross11, cross11], [cross12, cross12], alpha=0.3, color='royalblue', label=r"model1 $1\sigma$ C.I.")

    #plt.plot(lray2, chi2-209.417, "o", ms=4)
    plt.plot([60.201, 60.201], [-2, 10], "--", color='lightseagreen' ,label='model2 best-fit')
    plt.plot(dx2, dy2, '-')
    plt.fill_betweenx([0, 8], [cross21, cross21], [cross22, cross22], alpha=0.3, color='green', label=r'model2 $1\sigma$ C.I.')

    #plt.hlines(1, 50, 70, color="peru")

    plt.legend(loc='upper left')

    plt.xlabel(r"$L_{Ray}$/cm")
    plt.ylabel(r"$\Delta\chi^2$")
    plt.xlim(50, 70)
    plt.ylim(0, 8)
    plt.savefig("lray_profile.pdf")
    plt.show()
    """
