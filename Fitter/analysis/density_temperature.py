#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

def linearfunc(x, a, b):
    return a*x + b

T = [83.81, 86, 88, 90]
D =  [0.03549, 0.03513, 0.03481, 0.03449]
D1 = [0.03546, 0.03513, 0.03481, 0.03451]

from scipy.optimize import curve_fit
popt, pcov = curve_fit(linearfunc, T, D)
dx = np.arange(83, 91, 0.1)
dy = linearfunc(dx, *popt)
plt.plot(dx, dy, "-", color='blue')

plt.plot(T, D, "o", color='orange', label="Sinnock 1969")

popt1, pcov1 = curve_fit(linearfunc, T, D1)
dx1 = np.arange(83, 91, 0.1)
dy1 = linearfunc(dx1, *popt1)
plt.plot(dx1, dy1, "-", color='green')
plt.plot(T, D1, "s", color='darkviolet', label="NIST")

T_fit = 86.023
D_fit = 0.03507
plt.plot(T_fit, D_fit, "*", ms=9, color="red", label="from fitter")

plt.xlabel("temperature/K")
plt.ylabel(r"density mole/$cm^3$")

plt.tight_layout()
plt.grid(True)
plt.legend(loc='lower left')

plt.text(87, 0.0354, "y = %.6f*x + %.4f" %(popt[0], popt[1]),fontsize=13, color='blue')
plt.text(87, 0.0351, "y = %.6f*x + %.4f" %(popt1[0], popt1[1]),fontsize=13, color='green')

plt.savefig("den-T.pdf")
plt.show()