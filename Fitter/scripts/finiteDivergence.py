import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn-muted')

def babizcmodel(a0, aUV, aIR, ratio, x):
    A = a0 + aUV*x**2/(x**2-0.1066**2) + aIR*x**2/(x**2-0.9083**2)
    A *= ratio
    return np.sqrt(1+3*A/(3-A))


def finiteDiv(n, l):
    x = 17.5
    y = 11.6
    z = 36.5
    return ((x+y+z)/(x+y/n+z))**2


def main():
    corr_arr, wl_arr = [], []
    for wl in range(120, 150, 1):
        wl = wl/1000 #unit um
        wl_arr.append(wl)
        rindex = babizcmodel(0.3347, 0.0994, 0.0080, 0.03507/0.03449, wl)
        corr_arr.append( finiteDiv(rindex, wl) )

    plt.plot(wl_arr, corr_arr, "-", color='red')
    
    plt.xlabel("wavelengtth/nm")
    plt.ylabel("correction factor")

    plt.grid(True)
    plt.savefig("finiteDiv.pdf")
    plt.show()

if __name__ == '__main__':
    main()