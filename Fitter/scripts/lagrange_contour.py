#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':

    chimin = 207.998
    delta, lray, chi2 = [], [], []
    with open("../lag2d.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if float(data[1]) < 1000:
                delta.append(float(data[0]))
                lray.append(float(data[1]))
                chi2.append(float(data[2]))

minchi2 = 206.785
delta = np.array(delta)
lray = np.array(lray)
chi2 = np.array(chi2) - minchi2

from ROOT import TH2D
hist = TH2D("hist", "", 100, 1.36, 1.37, 100, 40, 80)
chi2_arr = []
for i in range(len(delta)):
    #x = int( (delta[i]-0.935)/ (0.97-0.935) * 100 ) + 1
    x = int( (delta[i]-1.36)/ (0.0001)) + 1
    y = int( (lray[i]-40) / 0.4) + 1
    hist.SetBinContent(x, y, chi2[i])

hist.SaveAs("hist2D.root")
#cont, edges = hist.numpy()
#xedges = edges[0][0]
#yedges = edges[0][1]
#plt.imshow(cont, interpolation="nearest", origin="lower", \
#    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

#plt.show()