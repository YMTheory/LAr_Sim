#!/usr/bin/env python
# coding=utf-8

from __future__ import division
import os
import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('seaborn-deep')

np.random.seed(1234)
import pymc
import scipy.stats as stats

def nominal():
    n = 100
    h = 61
    a = 2
    b = 2

    #p = pymc.Beta('p', alpha=a, beta=b)
    p = pymc.TruncatedNormal('p', mu=0.3, tau=10, a=0, b=1)
    y = pymc.Binomial('y', n=n, p=p, value=h, observed=True)
    m = pymc.Model([p, y])

    mc = pymc.MCMC(m, )
    mc.sample(iter=11000, burn=10000)
    plt.hist(p.trace(), 15, histtype='step', density=True, label='post');
    """
    x = np.linspace(0, 1, 100)
    plt.plot(x, stats.beta.pdf(x, a, b), label='prior');
    """
    a, b = plt.xlim()
    x = np.linspace(0, 1, 100)
    a, b = (0 - 0.3) / 0.1, (1 - 0.3) / 0.1
    plt.plot(x, stats.truncnorm.pdf(x, a, b, 0.3, 0.1), label='prior');
    plt.legend(loc='best');

if __name__ == '__main__':
    nominal()
    plt.show()
