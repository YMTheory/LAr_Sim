import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import ticker, cm

from ROOT import TH2D, TGraph, TCanvas, TFile, TAxis
def load_rfile(filename):
    ff = TFile(filename, "read")
    cc = ff.Get("c1")
    best = cc.GetPrimitive("best")
    hist = cc.GetPrimitive("hist")
    return best, hist

def graph_data(graph):
    x, y = [], []
    for i in range(graph.GetN()):            
        x.append(graph.GetPointX(i))
        y.append(graph.GetPointY(i))
    return np.array(x), np.array(y)

"""
N = 100
x = np.linspace(-3.0, 3.0, N)
y = np.linspace(-2.0, 2.0, N)

X, Y = np.meshgrid(x, y)
Z1 = np.exp(-(X)**2 - (Y)**2)
Z2 = np.exp(-(X * 10)**2 - (Y * 10)**2)
z = Z1 + 50 * Z2
z[:5, :5] = -1
print(z)

z = ma.masked_where(z <=0, z)
"""

if __name__ == '__main__':

    best, hist = load_rfile('/Users/yumiao/Documents/Works/LAr_Sim/Fitter/delta_ratio_profile2d_new.root')
    bestx, besty = graph_data(best)

    x, y = [], []
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()
    for i in range(1000):
        x.append(xaxis.GetBinCenter(i))
        y.append(yaxis.GetBinCenter(i))
    z = []
    for i in range(1000):
        tmp = []
        for j in range(1000):
            tmp.append(hist.GetBinContent(j, i))
        z.append(tmp)

    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots(figsize=(8, 6))
    cs = ax.contourf(X, Y, z, levels=[0, 2.3, 11.83])
    #cs = ax.contourf(X, Y, z, cmap=cm.PuBu_r, levels=[0, 2.3, 11.83])
    
    plt.plot([ 0.297, 0.317], [besty, besty], "-", color='orange')
    plt.plot([bestx, bestx], [0.932, 0.943], "-", color='orange')
    plt.plot(bestx, besty, "*", ms=8, color='hotpink', label="best fit(Model 1)")
    #plt.plot([ 126.505, 126.523 ], [besty, besty], "-", color='hotpink')
    #plt.plot([bestx, bestx], [ 140.104, 140.139 ], "-", color='hotpink')

    # Alternatively, you can manually set the levels
    # and the norm:
    # lev_exp = np.arange(np.floor(np.log10(z.min())-1),
    #                    np.ceil(np.log10(z.max())+1))
    # levs = np.power(10, lev_exp)
    # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())

    cbar = fig.colorbar(cs, orientation="vertical", pad=0.2)
    #cbar.ax.get_yaxis().set_ticks([])
    #cbar.ax.get_axis().labelpad = 15
    cbar.ax.set_yticks([0, 2.3, 11.83])
    cbar.ax.set_yticklabels([0, r"$1\sigma$", r"$3\sigma$"])
    #cbar.ax.set_ylabel('# of contacts', rotation=270)


    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    plt.grid(True)
    plt.xlim(0.24, 0.36)
    plt.ylim(0.90, 0.97)
    plt.xlabel(r"$\delta$")
    plt.ylabel('Xe absorption ratio')
    #plt.xlim(126.46, 126.57)
    #plt.ylim(140.02, 140.22)
    #plt.xlabel("absorption peak1 /nm")
    #plt.ylabel("absorption peak2 /nm")
    plt.legend()
    plt.grid(True)
    ax.set_axisbelow(True)
    plt.tight_layout()


    plt.savefig("profile2d_delta_ratio_new.pdf")
    plt.show()