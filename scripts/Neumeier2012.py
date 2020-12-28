#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import math

def loadData():
    # Draw Transmission from 2012
    wavelength = []
    trans = []
    trans_err = []
    with open('./G140ppb.txt') as f:
        for lines in f.readlines():
            line = lines.strip('\n')
            data = line.split(' ')
            if (float(data[0]) < 140):
                wavelength.append(float(data[0]))
                trans.append(float(data[1]))
                trans_err.append(float(data[2]))

    return wavelength, trans, trans_err

def calcAbs(wl, corr, trans, trans_err, L_Ray, L_Ray_err):
    T_att = trans * corr  # Fresnel correction
    A_att = -math.log(T_att)/math.log(10)
    A_att_err = 1/math.log(10)/T_att * corr * trans_err

    d = 5.8 #cm
    T_Ray = math.exp(-d/L_Ray)
    A_Ray = -math.log(T_Ray)/math.log(10)

    A_abs = A_att - A_Ray

    return A_abs, A_att_err

def fit1(xarr, yarr, yearr):
    # fitting with y = 1-exp(-a(x-b)) from literature
    # === fitting transmission spectrum === #
    from ROOT import TF1, TGraphErrors, TMath
    func = TF1("func", "1-TMath::Exp(-[0]*(x-[1]))", 118, 140)
    func.SetParameter(0, 0.234)
    func.SetParameter(1, 113.02)
    ge = TGraphErrors()
    for idd, x, y, yerr in zip([i for i in range(len(xarr))], xarr, yarr, yearr):
        ge.SetPoint(idd, x, y)
        ge.SetPointError(idd, 0, yerr)
    ge.Fit(func, "RE")
    drawx = []; drawy = []
    for i in range(118, 140, 1):
        drawx.append(i)
        drawy.append(func.Eval(i))
    return drawx, drawy

def rindex_zhou(l, a, b, c):
    # Zhou's Model
    A = 1.2055e-2*2/3.;
    rho_ratio = 34.49/(44.66e-3);
    l1 = 91.012;
    l2 = 89.892;
    l3 = 214.02;
    rindex = np.sqrt((3/(1-(A*rho_ratio*(a/(l1-1/l/l)+b/(l2-1/l/l)+c/(l3-1/l/l)))))-2)
    
    return rindex

def rindex_babizc(l, a0, aUV, aIR):
    # Babizc 2020
    lUV = 0.1066
    lIR = 0.9083

    A = a0 + aUV*l*l/(l**2-lUV**2) + aIR*l**2/(l**2-lIR**2)
    return np.sqrt(1+3*A/(3-A))

def raycalc(l, rindex):
    # transmission calculated from Rayleigh scattering length
    A = 1.2055e-2*2/3.;
    rho_ratio = 34.49/(44.66e-3);
    kT = 2.24442E-9
    kB = 1.380649E-23
    T = 90 # K
    f = 1E22

    delta = 0.3

    rayL = 1/ (8*(np.pi)**3/3/l**4 * ((rindex**2-1)*(rindex**2+2)/3)**2 * kB * kT *T *f * (6+3*delta)/(6-7*delta))
    #rayL = 1/ (8*(np.pi)**3/3/l**4 * ((rindex**2-1)*(rindex**2+2)/3)**2 * kB * kT *T *f)
    return rayL

def FresnelCorr(l, rindex):
    l = l*1000. #nm
    n_vac = 1; # vacuum
    a = 38 #eV2
    E0 = 12.8234 #eV
    gamma = 0.42357 #eV
    b = 247.188 #eV2
    E1 = 18.8448 #eV

    n_MgF2 = np.sqrt(1+a*(E0**2-(1240./l)**2)/((E0**2-(1240/l)**2)**2 + gamma**2*(1240./l)**2) + b/(E1**2-(1240/l)**2))

    factor = ( (1-((n_vac-n_MgF2)/(n_vac+n_MgF2))**2)/(1-((n_MgF2-rindex)/(n_MgF2+rindex))**2) )**2

    return factor


def main():

    plt.figure(0)

    wavelength, trans, trans_err = loadData()
    #plt.errorbar(wavelength, trans, yerr=trans_err, fmt='o', ms=2.5, lw=1.3, color='darkviolet', label='data')

    #drawx, drawy = fit1(wavelength, trans, trans_err)
    #plt.plot(drawx, drawy, "-", lw=1.6, color='forestgreen', label="fitting")

    rayx, rayT1, rayT2 = [], [], []
    corr_arr = []
    d = 5.8 # cm
    #p1, p2, p3 = 3.45227e-01, 1.89568e-02, 4.01598e+00  # Zhou's model fitting
    p1, p2, p3 = 3.35388e-01, 9.86866e-02, 8.04052e-03

    """
    for i in range(118, 140, 1):
        rayx.append(i)
        #rindex = rindex_zhou(i/1000., p1, p2, p3)
        rindex = rindex_babizc(i/1000., p1, p2, p3)
        corr = FresnelCorr(i/1000., rindex)
        corr_arr.append(corr)
        #rayT1.append(np.exp(-d/raycalc(i/1000., rindex))/corr)
        rayL = raycalc(i/1000., rindex)

    #plt.plot(rayx, rayT1, "--", lw=1.6, color='blue', label='Rayleigh: w/ Fresnel corr')
    """
    A_abs_arr, A_abs_err_arr = [], []
    for idx, td in enumerate(trans):
        wl = wavelength[idx]
        te = trans_err[idx]
        rindex = rindex_babizc(wl, p1, p2, p3)
        corr = FresnelCorr(wl, rindex)
        rayL = raycalc(wl, rindex)

        A_abs, A_abs_err = calcAbs(wl, corr, td, te, rayL, 0)
        A_abs_arr.append(A_abs)
        A_abs_err_arr.append(A_abs_err)

    plt.errorbar(wavelength, A_abs_arr, yerr=A_abs_err_arr,
                 fmt='o', ms=2.5, lw=1.3, color='darkviolet')

    plt.xlabel('wavalength/nm')
    #plt.ylabel('transmission')
    plt.ylabel('absorbance')
    plt.legend()
    plt.hlines(0, 118, 140, linestyles='dashed', color='blue')
    #plt.savefig('absorbance.pdf')

    #plt.figure(1)
    #plt.plot(rayx, corr_arr, "o-")
    #plt.title("Fresnel Correction")
    #plt.xlabel("wavelength/nm")
    #plt.ylabel("factor")
    #plt.grid(True)
    #plt.savefig('Fresnel90K.pdf')

    plt.show()


if __name__ == '__main__':
    main()
