import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import ticker, cm

from ROOT import TH2D

def read(filename):
    delta, ratio, chi = [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            delta.append(float(data[0]))
            ratio.append(float(data[1]))
            chi.append(float(data[2]))

    delta = np.array(delta)
    ratio = np.array(ratio)
    chi = np.array(chi)
    return delta, ratio, chi


def fill_th2d(delta, ratio, chi):
    hh = TH2D("profile2d", "", 500, 0, 0.5, 600, 0.88, 1.0)
    for i in range(len(delta)):
        binx = int(delta[i]/(0.5/500)) + 1
        biny = int((ratio[i]-0.88)/(0.12/600)) + 1
        print("%d, %d, %.3f" %(binx, biny, chi[i]))
        hh.SetBinContent(binx, biny, chi[i])

    return hh


def draw_cont(hist):
    x, y = [], []
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()
    for i in range(500):
        x.append(xaxis.GetBinCenter(i+1))
    for i in range(600):
        y.append(yaxis.GetBinCenter(i+1))
    z = []
    for i in range(600):
        tmp = []
        for j in range(500):
            #print(hist.GetBinContent(i+1, j+1))
            tmp.append(hist.GetBinContent(i+1, j+1))
        z.append(tmp)

    fig, ax = plt.subplots(figsize=(8, 6))
    #cs = ax.contourf(x, y, z, levels=[0, 2.3, 11.83])
    #cs = ax.contourf(X, Y, z, cmap=cm.PuBu_r, levels=[0, 2.3, 11.83])
    
    ax.imshow(z, extent=[x[0], x[-1], y[0], y[-1]])

    #cbar = fig.colorbar(cs, orientation='vertical', pad=0.2)
    plt.show()

def main():
    delta, ratio, chi = read("/scratchfs/higgs/yumiao/LAr/profile2D.txt")
    hh = fill_th2d(delta, ratio, chi)
    draw_cont(hh)

if __name__=="__main__":
    main()