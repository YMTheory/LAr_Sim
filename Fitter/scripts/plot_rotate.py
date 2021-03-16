#!/usr/bin/env python
# coding=utf-8

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

fig, ax =plt.subplot() 

# for rotate the axes and update.
for angle in range(0,360): 
    ax.view_init(30,angle)

plt.show()