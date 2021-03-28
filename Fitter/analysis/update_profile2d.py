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
            delta.append(float(data[0])*1000)
            ratio.append(float(data[1])*10000)
            chi.append(float(data[2]))

    delta = np.array(delta)
    ratio = np.array(ratio)
    chi = np.array(chi)
    chi = chi - 206.785
    return delta, ratio, chi


def fill_th2d(delta, ratio, chi):
    hh = TH2D("profile2d", "", 500, 0, 0.5, 600, 0.88, 1.0)
    for i in range(len(delta)):
        binx = int(delta[i]) + 1
        biny = int(((ratio[i]-8800)/2)) + 1
        #print("%.4f, %.4f, %d, %d, %.3f" %(delta[i], ratio[i], binx, biny, chi[i]))
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
            tmp.append(hist.GetBinContent(j+1, i+1))
        z.append(tmp)

    fig, ax = plt.subplots()
    ##cs = ax.contourf(x, y, z, levels=[0, 2.3, 11.83])
    cs = ax.contourf(x, y, z, levels=[0, 2.3, 11.83])
    #cs = ax.contourf(x, y, z, cmap=cm.PuBu_r, levels=[0, 2.3, 11.83])
    
    #cs = ax.imshow(z, extent=[x[0], x[-1], y[0], y[-1]], aspect='auto')
    cbar = fig.colorbar(cs, orientation='vertical', pad=0.15)
    cbar.ax.set_yticks([0, 2.3, 11.83])
    cbar.ax.set_yticklabels([0, r"$1\sigma$", r"$3\sigma$"])

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_axisbelow(True)

def main():
    filelist = []
    path = '../'
    for i in range(50):
        filename = "profile2D_delta+ratio_" + str(i) + ".txt"
        filelist.append(path+filename)
    delta, ratio, chi = [], [], []
    for ff in filelist:
        tmpdelta, tmpratio, tmpchi = read(ff)
        delta.extend(tmpdelta)
        ratio.extend(tmpratio)
        chi.extend(tmpchi)
    hh = fill_th2d(delta, ratio, chi)
    draw_cont(hh)

    plt.plot([0.2489, 0.3601], [0.937, 0.937], color='orange', lw=1)
    plt.plot([0.307, 0.307], [0.9218, 0.9533], color='orange', lw=1)
    plt.plot(0.307, 0.937, "*", color='hotpink', ms=8, label="Best Fit")

    plt.legend()
    plt.xlim(0, 0.5)
    plt.grid(True)
    plt.xlabel(r"$\delta$")
    plt.ylabel("Xe absorption peak ratio")
    plt.tight_layout()
    plt.savefig("profile2d_delta+ratio.pdf")
    plt.show()

if __name__=="__main__":
    main()