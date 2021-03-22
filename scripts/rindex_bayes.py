#!/usr/bin/env python
# coding=utf-8

import numpy as np

np.random.seed(1234)
import pymc

# observed data ...
xarr = np.array([0.3612, 0.3650, 0.4063, 0.4358, 0.4753, 0.5086, 0.5461, 0.5780, 0.6439])
yarr = np.array([1.2326, 1.2331, 1.2308, 1.2297, 1.2285, 1.2277, 1.2269, 1.2264, 1.2256])
yerr = np.array([0.001*i for i in yarr])

import matplotlib.pyplot as plt

plt.style.use('seaborn-muted')

#fig = plt.figure(figsize=(6, 6))
#plt.errorbar(xarr, yarr, yerr=yerr, fmt="o-", ms=5)
#plt.show()

import random
#transition_model = lambda x: [ x[0], random.gauss(x[1], 0.001), x[2] ]
transition_model = lambda x: [ random.gauss(x[0], 0.001), random.gauss(x[1], 0.001), random.gauss(x[2], 0.0001) ]
def transition_func(x):
    x1, x2 ,x0 = 0, 0, 0
    while x0<0.335*0.5 or x0>0.335*1.5:
        x0 = random.gauss(x[0], 0.005)
    while x1<0.099*0.5 or x1>0.099*1.5:
        x1 = random.gauss(x[1], 0.001)
    while x2<0.008*0.5 or x2>0.008*1.5:
        x2 = random.gauss(x[2], 0.0001)
    return [x0, x1, x2]
    

def prior(x) :
    m0, m1, m2 = 0.335, 0.099, 0.008
    s0, s1, s2 = 0.003, 0.003, 0.003
    return -(x[0]-m0)**2/2/s1**2 - (x[1]-m1)**2/2/s2**2 - (x[2]-m2)**2/2/s2**2  # gaussian prior
    #return 1.0

import math
def manual_log_like_normal(x, y, yerr, par):
    lik = 0
    for i in range(len(x)):
        A = par[0] + par[1] * x[i]**2/(x[i]**2 - 0.1066**2)  \
        + par[2]*x[i]**2/(x[i]**2 - 0.9083**2)
        if (1+3*A/(3-A)<0):
            pred = 10000000
        else:
            pred = math.sqrt(1+3*A/(3-A))
        lik += -(pred-y[i])**2/2/yerr[i]**2
    return lik

def acceptance(x, x_new):
    if x_new > x:
        return True
    else:
        accept = np.random.uniform(0, 1)
        return (accept < np.exp(x_new-x))

def metropolis_hastings(likelihood_computer,prior, transition_model, param_init,iterations,obsx, obsy, obse, acceptance_rule):
    p = param_init
    samples = []
    accepted = 0
    for i in range(iterations):
        p_new = transition_model(p)
        #p_new = transition_model(p)
        p_lik = likelihood_computer(obsx, obsy, obse, p)
        p_new_lik = likelihood_computer(obsx, obsy, obse, p_new)
        if (acceptance_rule(p_lik+prior(p), p_new_lik+prior(p_new))):
        #if (acceptance_rule(p_lik+np.log(prior(p)), p_new_lik+np.log(prior(p_new)))):
            p = p_new
            accepted += 1
        samples.append(p)
    print("Acceptance Ratio: %.4f" %(float(accepted)/float(iterations)) )
    return np.array(samples)


if __name__=='__main__':

    samples = metropolis_hastings(manual_log_like_normal, prior, transition_model, [0.4, 0.12, 0.012], 80000, xarr, yarr, yerr, acceptance)

    ax1 = plt.subplot(331)
    ax1.plot(samples[...,0], "-", lw=1)
    ax2 = plt.subplot(332)
    ax2.plot(samples[...,1], "-", lw=1)
    ax3 = plt.subplot(333)
    ax3.plot(samples[...,2], "-", lw=1)

    ax4 = plt.subplot(334)
    ax4.hist(samples[...,0], bins=40, histtype='step')
    ax5 = plt.subplot(335)
    ax5.hist(samples[...,1], bins=40, histtype='step')
    ax6 = plt.subplot(336)
    ax6.hist(samples[...,2], bins=40, histtype='step')

    ax7 = plt.subplot(337)
    ax7.scatter(samples[...,0], samples[...,1], s=0.1, color='peru')
    ax8 = plt.subplot(338)
    ax8.scatter(samples[...,0], samples[...,2], s=0.1, color='peru')
    ax9 = plt.subplot(339)
    ax9.scatter(samples[...,1], samples[...,2], s=0.1, color='peru')

    print("Parameter Estimated: %.4f+-%.4f, %.4f+-%.4f, %.4f+-%.4f" %(samples[...,0].mean(), np.std(samples[...,0]), \
        samples[...,1].mean(), np.std(samples[...,1]), samples[...,2].mean(), np.std(samples[...,2]) ) )

    plt.show()
