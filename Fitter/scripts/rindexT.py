#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

plt.style.use("seaborn-deep")

def rindex_ba(wl, s):
    l = wl
    lUV = 0.1066
    lIR = 0.9083
    a0 = 0.335
    aUV = 0.099
    aIR = 0.008
    A = s*( a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) )
    n = np.sqrt(1+3*A/(3-A))
    return n

m_wavelength = [0.3612, 0.3650, 0.4063, 0.4358, 0.4753, 0.5086, 0.5461, 0.5780, 0.6439] 
m_rindex0 = [1.2395, 1.2392, 1.2372, 1.2361, 1.2349, 1.2341, 1.2334, 1.2328, 1.2321] 
m_rindex1 = [1.2370, 1.2367, 1.2347, 1.2336, 1.2324, 1.2316, 1.2308, 1.2303, 1.2296]
m_rindex2 = [1.2349, 1.2346, 1.2326, 1.2315, 1.2303, 1.2295, 1.2287, 1.2282, 1.2274] 
m_rindex3 = [1.2326, 1.2331, 1.2308, 1.2297, 1.2285, 1.2277, 1.2269, 1.2264, 1.2256]
rerr = [0.001 for i in range(9)]

#plt.errorbar(m_wavelength, m_rindex0, yerr=rerr, fmt="o-", label="83.81K")
#plt.errorbar(m_wavelength, m_rindex1, yerr=rerr, fmt="o-", label="86K")
#plt.errorbar(m_wavelength, m_rindex2, yerr=rerr, fmt="o-", label="88K")
#plt.errorbar(m_wavelength, m_rindex3, yerr=rerr, fmt="o-", label="90K")

#dx = np.arange(0.3, 0.7, 0.001)
#dy3 = rindex_ba(dx, 1)
#dy0 = rindex_ba(dx, 3.549/3.449)
#dy1 = rindex_ba(dx, 3.513/3.449)
#dy2 = rindex_ba(dx, 3.481/3.448)
#plt.plot(dx, dy0, "--", color='darkviolet')
#plt.plot(dx, dy1, "--", color='darkviolet')
#plt.plot(dx, dy2, "--", color='darkviolet')
#plt.plot(dx, dy3, "--", color='darkviolet')

m_rindex0 = np.array(m_rindex0)
m_rindex1 = np.array(m_rindex1)
m_rindex2 = np.array(m_rindex2)
matrix = np.vstack([m_rindex0, m_rindex1, m_rindex2])
m_rindex = np.mean(matrix, axis=0)

m_diff0 = m_rindex0 - m_rindex
m_diff1 = m_rindex - m_rindex2
diff = np.vstack([m_diff0, m_diff1])
m_diff = np.max(diff, axis=0)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(m_wavelength, m_rindex, "o-")
ax.fill_between(m_wavelength, m_rindex-m_diff, m_rindex+m_diff, alpha=0.7)

plt.grid(True)
plt.xlabel("wavelength/um")
plt.ylabel("refindex")
#plt.legend()
plt.tight_layout()
plt.savefig("rindex_zone.pdf")

plt.show()
