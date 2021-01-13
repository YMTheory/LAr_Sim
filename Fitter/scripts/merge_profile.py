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

if __name__=='__main__':
    plt.style.use('seaborn-deep')

    graph1, zone1, line1 = load_graph('../delta_Ba.root')
    graph2, zone2, line2 = load_graph('../delta_Our.root')

    scanx1, scany1 = [], []
    for i in range(graph1.GetN()):
        scanx1.append(graph1.GetPointX(i))
        scany1.append(graph1.GetPointY(i))

    zonex1, zoney1 = [], []
    for i in range(zone1.GetN()):
        zonex1.append(zone1.GetPointX(i))
        zoney1.append(zone1.GetPointY(i))

    linex1, liney1 = [], []
    for i in range(line1.GetN()):
        linex1.append(line1.GetPointX(i))
        liney1.append(line1.GetPointY(i))

    plt.plot(scanx1, scany1, "-", color='royalblue')
    plt.plot(linex1, liney1, "--", color='blue', label="best fit(Model 1)")
    plt.fill_between([zonex1[0], zonex1[-1]], [0, 0], [zoney1[1], zoney1[1]], alpha=0.2, color='royalblue', label="5sigma C.I.(Model 1)")
    plt.ylim(0, liney1[-1])



    scanx2, scany2 = [], []
    for i in range(graph2.GetN()):
        scanx2.append(graph2.GetPointX(i))
        scany2.append(graph2.GetPointY(i))

    zonex2, zoney2 = [], []
    for i in range(zone1.GetN()):
        zonex2.append(zone2.GetPointX(i))
        zoney2.append(zone2.GetPointY(i))

    linex2, liney2 = [], []
    for i in range(line1.GetN()):
        linex2.append(line2.GetPointX(i))
        liney2.append(line2.GetPointY(i))

    plt.plot(scanx2, scany2, "-", color='forestgreen')
    plt.plot(linex2, liney2, "--", color='green', label="best fit (Model 2)")
    plt.fill_between([zonex2[0], zonex2[-1]], [0, 0], [zoney2[1], zoney2[1]], alpha=0.2, color='forestgreen', label="5sigma C.I.(Model 2)")

    plt.legend()
    plt.xlim(0, 0.4)
    plt.ylim(0, 400)

    plt.xlabel(r"$\delta$")
    plt.ylabel(r"$\Delta \chi^2$")

    plt.show()
