#!/bin/bash

import matplotlib.pyplot as plt
import numpy as np

def rindex_func(x):  # nm unit
    a = 38
    b = 247.188
    E0 = 12.8234
    E1 = 18.8448
    y = 0.42357
    E = 1240/x
    return np.sqrt( 1+ a*(E0**2-E**2)/((E0**2-E**2)**2+y**2*E**2) \
           + b/(E1**2-E**2) )

from ROOT import TGraph, TF1

if __name__ == "__main__":

    wavelength = [103.5, 104.65, 105.25, 106.00, 106.95, 108.15, 110.25, 111.77, 114.55, 117.85, 121.53, 130.15, 140.00, 150.00, 170.00, 200.00]
    expdata    = [1.9710, 1.9140, 1.8875, 1.8600, 1.8325, 1.7980, 1.7505, 1.7245, 1.6840, 1.6490, 1.6160, 1.5650, 1.5255, 1.4985, 1.4705, 1.4450]
    wavelength = np.array(wavelength)

    plt.plot(wavelength, expdata, 'o-', ms=5, label="LiF 193deg exp")
    plt.plot(wavelength, rindex_func(wavelength), label="Parameterization")

    plt.xlabel("wavelength/nm")
    plt.ylabel("rindex")
    plt.legend()
    plt.show()