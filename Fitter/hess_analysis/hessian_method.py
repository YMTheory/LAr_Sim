#!/usr/bin/env python
# coding=utf-8

import numpy as np

best_a0 = 0.330269
best_aUV = 0.103343
best_aIR = 0.0076295
best_delta = 0.249214

def load_cov(filename):
    with open(filename) as f:
        cov = []
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            tmp = []
            for elem in data:
                tmp.append(float(elem))
            cov.append(tmp)

    cov = np.array(cov)
    return cov

def rindex_ba(wl, a0, aUV, aIR):
    l = wl
    lUV = 0.1066
    lIR = 0.9083

    A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ;
    n = np.sqrt(1+3*A/(3-A))
    return n

def finite_diff_rindex(idx, step):
    a0 = best_a0
    aUV = best_aUV
    aIR = best_aIR
    rindex = rindex_ba(0.128, a0, aUV, aIR)
    if idx == 0: # a0
        delta_rindex1 = rindex_ba(0.128, a0-step, aUV, aIR)-rindex
        delta_rindex2 = rindex_ba(0.128, a0+step, aUV, aIR)-rindex
    if idx == 1: # aUV
        delta_rindex1 = rindex_ba(0.128, a0, aUV-step, aIR)-rindex
        delta_rindex2 = rindex_ba(0.128, a0, aUV+step, aIR)-rindex
    if idx == 2: # aIR
        delta_rindex1 = rindex_ba(0.128, a0, aUV, aIR-step)-rindex
        delta_rindex2 = rindex_ba(0.128, a0, aUV, aIR+step)-rindex
        
    print("%.7f, %.7f, %.7f" %(delta_rindex1, delta_rindex2, step))
    dev = (delta_rindex2 - delta_rindex1) / 2 /step

    return dev


def raylength(wl, delta, n):
    l = wl
    kT = 2.24442E-9
    kB = 1.380649E-23
    T = 90
    f = 1e22;
    rayL = 1 / (8*np.pi**3/3/l**4 * ((n**2-1)*(n**2+2)/3)**2 * kT * kB * T * f *(6+3*delta)/(6-7*delta))

    return rayL


def finite_diff_rayL(idx, step):
    a0 = best_a0
    aUV = best_aUV
    aIR = best_aIR
    delta = best_delta
    rayL = raylength(0.128, delta, rindex_ba(0.128, a0, aUV, aIR))
    if idx == 0: # a0
        rayL1 = raylength(0.128, delta, rindex_ba(0.128, a0-step, aUV, aIR)) - rayL
        rayL2 = raylength(0.128, delta, rindex_ba(0.128, a0+step, aUV, aIR)) - rayL
    if idx == 1: # aUV
        rayL1 = raylength(0.128, delta, rindex_ba(0.128, a0, aUV-step, aIR)) - rayL
        rayL2 = raylength(0.128, delta, rindex_ba(0.128, a0, aUV+step, aIR)) - rayL
    if idx == 2: # aIR
        rayL1 = raylength(0.128, delta, rindex_ba(0.128, a0, aUV, aIR-step)) - rayL
        rayL2 = raylength(0.128, delta, rindex_ba(0.128, a0, aUV, aIR+step)) - rayL
    if idx == 3: # delta
        rayL1 = raylength(0.128, delta-step, rindex_ba(0.128, a0, aUV, aIR)) - rayL
        rayL2 = raylength(0.128, delta+step, rindex_ba(0.128, a0, aUV, aIR)) - rayL

        
    dev = (rayL2 - rayL1) / 2 /step
    return dev



import math
def main():
    cov = load_cov("./cov_matrix.txt")

    dev = []
    dev.append(finite_diff_rayL(0, 0.0001))
    dev.append(finite_diff_rayL(1, 0.0001))
    dev.append(finite_diff_rayL(2, 0.0001))
    dev.append(finite_diff_rayL(3, 0.0001))
    print(dev)

    delta_chi2 = 1
    delta_X2 = 0
    for i in range(4):
        for j in range(4):
            delta_X2 += delta_chi2 * (dev[i]*dev[j]*cov[i, j])

    print("final deviation for rindex@128nm : %.6f" %math.sqrt(delta_X2) )


from ROOT import TFile, TCanvas, TGraph
def load_graph(filename):
    ff = TFile(filename, "read")
    cc = ff.Get("c1")
    graph = cc.GetPrimitive("graph")
    zone = cc.GetPrimitive("zone")
    line = cc.GetPrimitive("line")
    return graph, zone, line

def approx_check():
    cov = load_cov("./hessian_matrix.txt")

    graph, zone, line = load_graph("../outputs/aIR_model1.root")
    scanx1, scany1, scanl1 = [], [], []
    for i in range(graph.GetN()):
        scanx1.append(graph.GetPointX(i))
        scany1.append(graph.GetPointY(i))

    plt.plot(scanx1, scany1, "-", color='royalblue', label="calculation")

    hess_pred = []
    for i in scanx1:
        hess_pred.append(cov[2, 2]*(i-best_aIR) * (i-best_aIR))

    plt.plot(scanx1, hess_pred, "--", color='lightseagreen', label="Hessian approx")

    #plt.xlim(0.2, 0.3)
    #plt.xlim(0.099, 0.108)
    plt.xlim(0.002, 0.014)
    plt.ylim(0, 30)
    
    plt.legend()
    plt.xlabel("aIR")
    #plt.xlabel(r"$\delta$")
    plt.ylabel(r"$\Delta \chi^2$")
    plt.grid(True)
    plt.savefig("aIR_approx.pdf")
    plt.show()

import matplotlib.pyplot as plt

if __name__ == '__main__':
    plt.style.use('seaborn-deep')
    #main()

    approx_check()
            
