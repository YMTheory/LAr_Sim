import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def read_txt(filename):
    lr, chi2origin, n128 = [], [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if data[0] == "factor":
                lr.append(float(data[1]))
            else:
                n128.append(float(data[0]))
                chi2origin.append(float(data[1]) - lr[-1]*n128[-1])

    lr = np.array(lr)
    n128 = np.array(n128)
    chi2origin = np.array(chi2origin)

    return lr, n128, chi2origin


def func(x, a, b, c):
    return a*x**2 + b*x + c


def calculate_interval(val, a, b, c):
    x1 = (-b-np.sqrt(b**2-4*a*(c-val))) / 2/a
    x2 = (-b+np.sqrt(b**2-4*a*(c-val))) / 2/a
    return x1, x2

def find_interval(x, y):
    val1, val2 = 0, 0
    for i in range(len(x)-1):
        if y[i] >=1 and y[i+1]<=1:
            val1 = np.interp(1, [y[i], y[i+1]], [x[i], x[i+1]])
        if y[i] <=1 and y[i+1]>=1:
            val2 = np.interp(1, [y[i], y[i+1]], [x[i], x[i+1]])
    return val1, val2


if __name__ == "__main__" :
    lr, n128, chi2origin = read_txt("lagrange_rindex128.txt")
    #lr, n128, chi2origin = read_txt("lagrange_lray128.txt")
    nn = len(n128)
    n128_sorted = n128
    chi2origin_sorted = chi2origin
    suffixStart = 0
    while suffixStart != nn:
        for i in range(suffixStart, nn):
            if n128_sorted[i] < n128_sorted[suffixStart]:
                n128_sorted[suffixStart] ,n128_sorted[i] = n128_sorted[i], n128_sorted[suffixStart]
                chi2origin_sorted[suffixStart] ,chi2origin_sorted[i] = chi2origin_sorted[i], chi2origin_sorted[suffixStart]
        suffixStart += 1
    
    xleft, xright = find_interval(n128_sorted, chi2origin_sorted-np.min(chi2origin_sorted))
    print(xleft, xright)

    popt, pcov = curve_fit(func, n128_sorted, chi2origin_sorted-np.min(chi2origin_sorted))
    xleft, xright = calculate_interval(1, *popt)
    print(xleft, xright)
    

    ## Fitting parabola
    ##dx = np.linspace(np.min(n128), np.max(n128), 100)
    ##dy = func(dx, *popt)
    ### part of the data
    #n128_draw = [n128_sorted[i] for i in range(0, len(n128_sorted)-5, 5)]
    #chi2origin_draw = [chi2origin_sorted[i] for i in range(0, len(chi2origin_sorted)-5, 5)]
    ## interval
    ##xleft, xright = calculate_interval(1, *popt)
    ##print(xleft, xright)
    #
    #fig, ax = plt.subplots()
    #ax.plot(n128_draw, chi2origin_draw- np.min(chi2origin), "o", ms=5, color="red", zorder=1)
    ##ax.plot(dx, dy, "-", color="black")
    #ax.hlines(1, np.min(n128), np.max(n128), color="black", lw=1.5)
    #ax.fill_betweenx([0, np.max(chi2origin-np.min(chi2origin))], [xleft, xleft], [xright, xright], color="blue", alpha=0.3)
    ##ax.set_xlabel("n(128) nm", fontsize=14)
    #ax.set_xlabel(r"L$_\mathrm{Ray}$(128 nm) [cm]", fontsize=14)
    #ax.set_ylabel(r"$\Delta \chi^2$", fontsize=14)
    #plt.tight_layout()
    ##plt.savefig("refractive_index_128nm_lagrange.pdf")
    #plt.savefig("lray_128nm_lagrange.pdf")
    #plt.show()



