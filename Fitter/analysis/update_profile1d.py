#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def read(filename):
    var, chi2 = [], []
    chimin = 10000
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            var.append(float(data[0]))
            chi2.append(float(data[1]))
            if chimin > chi2[-1]:
                chimin = chi2[-1]
    var = np.array(var)
    chi2 = np.array(chi2)
    return var, chi2, chimin

def find_interval(var, chi):
    left_margin, right_margin = 0, 0
    for i in range(len(var)-1):
        if chi[i] >= 1 and chi[i+1] <= 1:
            left_margin = var[i+1] - (var[i+1]-var[i]) * (1 - chi[i+1]) / (chi[i] - chi[i+1])
        if chi[i] <= 1 and chi[i+1] >= 1:
            right_margin = var[i] + (var[i+1]-var[i]) * (1 - chi[i]) / (chi[i+1] - chi[i])

    return left_margin, right_margin



def main():

    delta, chi, chimin = read("./profile1D_ratio.txt")
    left, right = find_interval(delta, chi-chimin)
    print("[%.4f, %.4f]" %(left, right))

    plt.plot(delta, chi-chimin, "-", color='blue')
    plt.fill_betweenx([0, 20], left, right, color='violet', alpha=0.4)
    
    #plt.xlabel(r"$\delta$")
    plt.ylabel(r"$\Delta\chi^2$")
    #plt.xlim(0, 0.5)

    plt.grid(True)
    #plt.savefig("profile1D_delta.pdf")
    plt.show()


if __name__ == "__main__":
    main()

