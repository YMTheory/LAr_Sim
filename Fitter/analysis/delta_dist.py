#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

delta_arr = []
with open("../delta_toyMC.txt") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        delta_arr.append(float(line))

plt.hist(delta_arr, bins=100, range=(0, 0.5), alpha=0.5)

plt.tight_layout()
plt.xlabel("delta")
plt.show()