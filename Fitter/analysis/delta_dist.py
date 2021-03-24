#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

delta_arr = []
a0_arr, aUV_arr, aIR_arr = [], [], []

with open("../apar_toyMC.txt") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        #delta_arr.append(float(line))
        a0_arr.append(float(data[0]))
        aUV_arr.append(float(data[1]))
        aIR_arr.append(float(data[2]))


#plt.hist(delta_arr, bins=100, range=(0.3335, 0.336), alpha=0.5, edgecolor="black")
#plt.hist(aUV_arr, bins=100, range=(0.098, 0.101), alpha=0.5, edgecolor="black")
plt.hist(aIR_arr, bins=100, range=(0.0025, 0.0125), alpha=0.5, edgecolor="black")

plt.tight_layout()
plt.xlabel(r"$a_{UV}$")
plt.savefig("aIR_dist.pdf")
plt.show()