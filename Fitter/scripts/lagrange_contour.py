#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':

    chimin = 207.998
    delta, lray, chi2 = [], [], []
    with open("../log2d_rayL+rindex.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            delta.append(float(data[0]))
            lray.append(float(data[1]))
            chi2.append(float(data[2]))

delta = np.array(delta)
lray = np.array(lray)
chi2 = np.array(chi2)

from ROOT import TH2D
hist = TH2D("hist", "", 100, 50, 130, 100, 1.3, 1.4)
chi2_arr = []
for i in range(len(delta)):
    #x = int( (delta[i]-0.935)/ (0.97-0.935) * 100 ) + 1
    x = int( (delta[i]-1.3)/ (0.1) * 100 ) + 1
    y = int( (lray[i]-50) / 80*100 ) + 1
    hist.SetBinContent(x, y, chi2[i])

hist.SaveAs("hist2D.root")
#cont, edges = hist.numpy()
#xedges = edges[0][0]
#yedges = edges[0][1]
#plt.imshow(cont, interpolation="nearest", origin="lower", \
#    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

#plt.show()