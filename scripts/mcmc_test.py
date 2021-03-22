#!/usr/bin/env python
# coding=utf-8

from __future__ import division
import os
import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

plt.style.use('seaborn-deep')

from mpl_toolkits.mplot3d import Axes3D
import scipy.stats as stats
from functools import partial

np.random.seed(1234)

def target(lik, prior, n, h, theta):
    if theta < 0 or theta > 1:
        return 0.0
    else:
        return lik(n, theta).pmf(h) * prior.pdf(theta)

def onecoin():
    n = 100
    h = 61
    a = 10
    b = 10
    lik = st.binom
    prior = st.beta(a, b)
    sigma = 0.3 
    thetas = np.linspace(0, 1, 200)

    naccept = 0
    theta = 0.1
    niters = 10000
    samples = np.zeros(niters+1)
    samples[0] = theta
    for i in range(niters):
        theta_p = theta + st.norm(0, sigma).rvs()
        rho = min(1, target(lik, prior, n, h, theta_p)/target(lik, prior, n, h, theta))
        u = np.random.uniform()
        if u < rho:
            naccept += 1
            theta = theta_p
        samples[i+1] = theta
    nmcmc = len(samples)//2
    print('Efficiency = %.4f' %(naccept/niters))

    post = st.beta(h+a, n-h+b)

    plt.figure(0)
    plt.hist(samples[nmcmc:],  40, histtype='step', density=True, linewidth=1, label="Distribution of prior samples")
    plt.hist(prior.rvs(nmcmc), 40, histtype='step', density=True, linewidth=1, label='Distribution of posterior samples')
    plt.plot(thetas, post.pdf(thetas), c='red', linestyle='--', alpha=0.5, label='True posterior')
    plt.figure(1)
    plt.plot(samples, "-")

    plt.legend()
    plt.show()

def mh_coin(niters, n, h, theta, lik, prior, sigma):
    samples = [theta]
    while len(samples) < niters:
        theta_p = theta + st.norm(0, sigma).rvs()
        rho = min(1, target(lik, prior, n, h, theta_p)/target(lik, prior, n, h, theta ))
        u = np.random.uniform()
        if u < rho:
            theta = theta_p
        samples.append(theta)
    return samples

def convergence():
    n = 100
    h = 61
    a = 10
    b = 10
    lik = st.binom
    prior = st.beta(a, b)
    sigma = 0.05
    niters = 100

    sampless = [mh_coin(niters, n, h, theta, lik, prior, sigma) for theta in np.arange(0.1, 1, 0.2)]

    for samples in sampless:
        plt.plot(samples, '-o')
    plt.xlim([0, niters])
    plt.ylim([0, 1])
    plt.show()


def bern(theta, z, N):
    return np.clip(theta**z*(1-theta)**(N-z), 0, 1)

def bern2(theta1, theta2, z1, z2, N1, N2):
    """Bernoulli likelihood with N trials and z successes."""
    return bern(theta1, z1, N1) * bern(theta2, z2, N2)

def make_thetas(xmin, xmax, n):
    xs = np.linspace(xmin, xmax, n)
    widths = (xs[1:] - xs[:-1]) / 2.0
    thetas = xs[:-1] + widths
    return thetas

def make_plots(X, Y, prior, likelihood, posterior, projection=None):
    fig, ax = plt.subplots(1,3, subplot_kw=dict(projection=projection), figsize=(12,3))
    #fig, ax = plt.subplots(1,3, subplot_kw=dict(projection=projection, aspect='equal'), figsize=(12,3))
    if projection == '3d':
        ax[0].plot_surface(X, Y, prior, alpha=0.3, cmap=plt.cm.jet)
        ax[1].plot_surface(X, Y, likelihood, alpha=0.3, cmap=plt.cm.jet)
        ax[2].plot_surface(X, Y, posterior, alpha=0.3, cmap=plt.cm.jet)
    else:
        ax[0].contour(X, Y, prior)
        ax[1].contour(X, Y, likelihood)
        ax[2].contour(X, Y, posterior)
    ax[0].set_title('Prior')
    ax[1].set_title('Likelihood')
    ax[2].set_title('Posteior')
    plt.tight_layout()

def twocoins():
    thetas1 = make_thetas(0, 1, 101)
    thetas2 = make_thetas(0, 1, 101)
    X, Y = np.meshgrid(thetas1, thetas2)
    thetas1 = make_thetas(0, 1, 101)

    a = 2
    b = 3
    z1 = 11
    N1 = 14
    z2 = 7
    N2 = 14

    prior = st.beta(a, b).pdf(X) * st.beta(a, b).pdf(Y)
    likelihood = bern2(X, Y, z1, z2, N1, N2)
    posterior = st.beta(a+z1, b+N1-z1).pdf(X) * st.beta(a+z2, b+N2-z2).pdf(Y)

    make_plots(X, Y, prior, likelihood, posterior)
    make_plots(X, Y, prior, likelihood, posterior, projection='3d')
    plt.show()

if __name__ == '__main__':

    onecoin()
    #convergence()