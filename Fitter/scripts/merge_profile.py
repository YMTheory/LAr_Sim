#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

from ROOT import TFile, TCanvas, TGraph
def load_graph(filename):
    ff = TFile(filename, "read")
    cc = ff.Get("c1")
    graph = cc.GetPrimitive("graph")
    zone = cc.GetPrimitive("zone")
    line = cc.GetPrimitive("line")
    return graph, zone, line

def rindex_ba(wl):
    l = wl
    lUV = 0.1066
    lIR = 0.9083
    a0 = 0.335
    aUV = 0.099
    aIR = 0.008
    A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ;
    n = np.sqrt(1+3*A/(3-A))
    return n

def rindex_our(wl):
    A = 1.2055e-2*2/3.;
    rho_ratio = 34.49/(44.66e-3);
    l1 = 91.012;
    l2 = 89.892;
    l3 = 214.02;
    l = wl
    a = 0.280
    b = 0.079
    c = 4.025

    return np.sqrt((3/(1-(A*rho_ratio*(a/(l1-1/l/l)+b/(l2-1/l/l)+c/(l3-1/l/l)))))-2);


def raylength(wl, delta, n):
    l = wl
    kT = 2.24442E-9
    kB = 1.380649E-23
    T = 90
    f = 1e22;
    rayL = 1 / (8*np.pi**3/3/l**4 * ((n**2-1)*(n**2+2)/3)**2 * kT * kB * T * f *(6+3*delta)/(6-7*delta))

    return rayL


if __name__=='__main__':
    plt.style.use('seaborn-deep')

    #graph1, zone1, line1 = load_graph('/Users/yumiao/Documents/Works/LAr_Sim/Fitter/delta.root')
    graph1, zone1, line1 = load_graph('/Users/yumiao/Documents/Works/LAr_Sim/Fitter/outputs/delta_5sigma.root')
    #graph2, zone2, line2 = load_graph('/Users/yumiao/Documents/Works/LAr_Sim/Fitter/outputs/delta_model2.root')

    scanx1, scany1, scanl1 = [], [], []
    for i in range(graph1.GetN()):
        scanx1.append(graph1.GetPointX(i))
        scany1.append(graph1.GetPointY(i))
        rindex = rindex_ba(0.128)
        scanl1.append(raylength(0.128, scanx1[-1], rindex))

    zonex1, zoney1, zonel1 = [], [], []
    for i in range(zone1.GetN()):
        zonex1.append(zone1.GetPointX(i))
        zoney1.append(zone1.GetPointY(i))
        rindex = rindex_ba(0.128)
        zonel1.append(raylength(0.128, zonex1[-1], rindex))

    linex1, liney1, linel1 = [], [], []
    for i in range(line1.GetN()):
        linex1.append(line1.GetPointX(i))
        liney1.append(line1.GetPointY(i))
        rindex = rindex_ba(0.128)
        linel1.append(raylength(0.128, linex1[-1], rindex))

    fig = plt.figure(figsize=(6, 4))
    figr = fig.add_subplot(111)

    plt.plot(scanx1, scany1, "-", color='royalblue')
    plt.plot(linex1, liney1, "--", color='blue', label="best fit(Model 1)")
    plt.fill_between([zonex1[0], zonex1[-1]], [0, 0], [zoney1[1], zoney1[1]], alpha=0.2, color='purple', label="5sigma C.I.(Model 1)")

    #plt.hlines(25, 0, 0.5, color='magenta', linestyles="-.")
    #plt.text(0.04, 60, r"$5 \sigma$ C.I.")
    #plt.hlines(25, 0.9, 1.0, color='magenta', linestyles="-.")
    #plt.text(0.91, 30, r"$5 \sigma$ C.I.")
    #plt.plot(scanl1, scany1, '-', color='royalblue')
    #plt.plot(linel1, liney1, "--", color='blue', label="best fit(Model 1)")
    #plt.fill_between([zonel1[0], zonel1[-1]], [0, 0], [zoney1[1], zoney1[1]], alpha=0.2, color='royalblue', label="5sigma C.I.(Model 1)")
    plt.ylim(0, liney1[-1])


    """
    scanx2, scany2, scanl2 = [], [], []
    for i in range(graph2.GetN()):
        scanx2.append(graph2.GetPointX(i))
        scany2.append(graph2.GetPointY(i))
        rindex = rindex_our(0.128)
        scanl2.append(raylength(0.128, scanx2[-1], rindex))

    zonex2, zoney2, zonel2 = [], [], []
    for i in range(zone1.GetN()):
        zonex2.append(zone2.GetPointX(i))
        zoney2.append(zone2.GetPointY(i))
        rindex = rindex_our(0.128)
        zonel2.append(raylength(0.128, zonex2[-1], rindex))

    linex2, liney2, linel2 = [], [], []
    for i in range(line1.GetN()):
        linex2.append(line2.GetPointX(i))
        liney2.append(line2.GetPointY(i))
        rindex = rindex_our(0.128)
        linel2.append(raylength(0.128, linex2[-1], rindex))

    #plt.plot(scanx2, scany2, "-", color='forestgreen')
    #plt.plot(linex2, liney2, "--", color='green', label="best fit (Model 2)")
    #plt.fill_between([zonex2[0], zonex2[-1]], [0, 0], [zoney2[1], zoney2[1]], alpha=0.2, color='forestgreen', label="5sigma C.I.(Model 2)")
    plt.plot(scanl2, scany2, "-", color='forestgreen')
    plt.plot(linel2, liney2, "--", color='green', label="best fit (Model 2)")
    plt.fill_between([zonel2[0], zonel2[-1]], [0, 0], [zoney2[1], zoney2[1]], alpha=0.2, color='forestgreen', label="5sigma C.I.(Model 2)")
    """

    plt.legend(loc='upper left')
    #plt.xlim(0.92, 0.96)
    plt.xlim(0.0, 0.50)
    #plt.ylim(0, 20)
    #plt.xlim(126.45, 126.57)
    #plt.xlim(140.07, 140.17)

    plt.xlabel(r"$\delta$")
    #plt.xlabel("Xe absorption peak ratio")
    plt.ylabel(r"$\Delta \chi^2$")
    plt.grid(True)
    #figr.axes.get_xaxis().set_visible(False)
    #plt.xticks([])

    plt.savefig("profile_delta_0zone_model1.pdf")
    plt.show()
