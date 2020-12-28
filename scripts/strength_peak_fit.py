#!/usr/bin/env python
# coding=utf-8

import matplotlib.pyplot as plt
import numpy as np

if __name__ =="__main__":
    
    T = [180.0, 163.5, 158.7, 155.3, 148.7, 140.1, 126.3, 111.8, 104.4, 86.8]
    rho = [4.17, 6.29, 8.16, 9.9, 13, 14.6, 16.9, 18.5, 19.5, 21.1]
    fa = [0.029, 0.037, 0.041, 0.067, 0.076, 0.069, 0.067, 0.074, 0.063, 0.054]
    Ea = [11.59, 11.58, 11.58, 11.58, 11.61, 11.63, 11.67, 11.71, 11.73, 11.79]
    fb = [0.206, 0.183, 0.194, 0.196, 0.181, 0.180, 0.188, 0.162, 0.167, 0.207]
    Eb = [11.85, 11.84, 11.86, 11.87, 11.89, 11.89, 11.94, 11.97, 11.98, 12.05]

    plt.subplot(221)
    plt.plot(T, fa, "o", label='fa')
    plt.legend()
    
    plt.subplot(222)
    plt.plot(T, fb, "o", label='fb')
    plt.xlabel('T/K')
    plt.legend()

    plt.subplot(223)
    plt.plot(T, fa, "o", label='Ea')
    plt.legend()

    plt.subplot(224)
    plt.plot(T, fa, "o", label='Eb')
    plt.xlabel('T/K')
    plt.legend()

    plt.show()
