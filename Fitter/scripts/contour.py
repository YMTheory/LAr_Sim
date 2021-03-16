#!/usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

from ROOT import TGraph, TCanvas, TFile
def load_graph(filename):
    ff = TFile(filename, "read")
    cc = ff.Get("c1")
    best = cc.GetPrimitive("best")
    sigma1 = cc.GetPrimitive("sigma1")
    sigma5 = cc.GetPrimitive("sigma5")
    return best, sigma1, sigma5

def graph_data(graph):
    x, y = [], []
    for i in range(graph.GetN()):            
        x.append(graph.GetPointX(i))
        y.append(graph.GetPointY(i))
    return np.array(x), np.array(y)


import numpy.linalg as linalg

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  linalg.eig(np.dot(linalg.inv(S), C))
    n =  np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def find_ellipse(x, y):
    xmean = x.mean()
    ymean = y.mean()
    x = x - xmean
    y = y - ymean
    a = fitEllipse(x,y)
    center = ellipse_center(a)
    center[0] += xmean
    center[1] += ymean
    phi = ellipse_angle_of_rotation(a)
    axes = ellipse_axis_length(a)
    x += xmean
    y += ymean
    return center, phi, axes

def rindex_ba(wl):
    l = wl
    lUV = 0.1066
    lIR = 0.9083
    a0 = 0.335
    aUV = 0.0987
    aIR = 0.008
    A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ;
    n = np.sqrt(1+3*A/(3-A))
    return n

def rindex_our(wl):
    A = 1.2055e-2*2/3.;
    rho_ratio = 34.49/(44.66e-3);
    l1 = 91.012;
    l2 = 89.892;
    l3 = 214.02;
    l = wl
    a = 2.21678e-04
    b = 3.38135e-01
    c = 4.06840

    return np.sqrt((3/(1-(A*rho_ratio*(a/(l1-1/l/l)+b/(l2-1/l/l)+c/(l3-1/l/l)))))-2);



def raylength(wl, delta, n):
    l = wl
    kT = 2.24442E-9
    kB = 1.380649E-23
    T = 90
    f = 1e22;
    rayL = 1 / (8*np.pi**3/3/l**4 * ((n**2-1)*(n**2+2)/3)**2 * kT * kB * T * f *(6+3*delta)/(6-7*delta))

    return rayL

if __name__=='__main__':
   
    
    best, sigma1, sigma5 = load_graph("/Users/yumiao/Documents/Works/LAr_Sim/Fitter/outputs/contour_model1.root")
    bestx, besty = graph_data(best)
    sigma1x, sigma1y = graph_data(sigma1)
    sigma5x, sigma5y = graph_data(sigma5)

    #best_our, sigma1_our, sigma5_our = load_graph("/Users/yumiao/Documents/Works/LAr_Sim/Fitter/outputs/contour_model2.root")
    #bestx_our, besty_our = graph_data(best_our)
    #sigma1x_our, sigma1y_our = graph_data(sigma1_our)
    #sigma5x_our, sigma5y_our = graph_data(sigma5_our)


    # Model 1
    points = []
    for x,y in zip(sigma1x, sigma1y):
        points.append([x, y])

    a_points = np.array(points)
    x1 = a_points[:, 0]
    y1 = a_points[:, 1]
    center1, phi1, axes1 = find_ellipse(x1, y1)
    print("axes length: %.3f, %.3f" %(axes1[0], axes1[1]))
    print("ratation angle: %.5f" %( float(phi1)/np.pi*180) )

    xx1, yy1, zz1 = [], [], []
    for i in range(1000):
        ang = np.pi*2/1000.*i
        xx1.append( axes1[0]*np.cos(ang)*np.cos(phi1)-axes1[1]*np.sin(ang)*np.sin(phi1)+center1[0]  )
        yy1.append( axes1[0]*np.cos(ang)*np.sin(phi1)+axes1[1]*np.sin(ang)*np.cos(phi1)+center1[1]  )
        zz1.append(raylength(0.128, yy1[-1], rindex_ba(0.128)))


    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(7.0, 6))

    plt.plot(bestx, besty, "*", ms=8, color='seagreen', label="best fit(Model 1)")
    plt.plot([0.2557, 0.2776], [besty, besty], "-", color='seagreen')
    plt.plot([bestx, bestx], [0.9322, 0.9432], "-", color='seagreen')
    #plt.plot(x1, y1, "o", ms=4, color='green')
    plt.plot(xx1, yy1, "-", color='lightseagreen', label='1sigma allowed region (Model 1)')
    #print("Model1 delta = %.3f, rindex@128nm = %.3f, raylength@128nm = %.3f cm" %(besty, rindex_ba(0.128), raylength(0.128, besty, rindex_ba(0.128))) )
    #bestz = raylength(0.128, besty, rindex_ba(0.128))
    #plt.plot(bestx, bestz, "*", ms=8, color='seagreen', label="best fit(Model 1)")
    #plt.plot(xx1, zz1, "-", color='lightseagreen', label='1sigma allowed region(Model 1)')

    #import turtle
    #turtle.shape("circle")
    #turtle.shapesize(0.024, 0.008)
    #turtle.setpos(bestx, besty)
    #turtle.fillcolor('pink')

    """
    # Model 2
    points = []
    for x,y in zip(sigma1x_our, sigma1y_our):
        points.append([x, y])

    a_points = np.array(points)
    x2 = a_points[:, 0]
    y2 = a_points[:, 1]
    center2, phi2, axes2 = find_ellipse(x2, y2)
    print("axes length: %.3f, %.3f" %(axes2[0], axes2[1]))
    print("ratation angle: %.5f" %( float(phi2)/np.pi*180) )

    xx2, yy2, zz2 = [], [], []
    for i in range(1000):
        ang = np.pi*2/1000.*i
        xx2.append( axes2[0]*np.cos(ang)*np.cos(phi2)-axes2[1]*np.sin(ang)*np.sin(phi2)+center2[0]  )
        yy2.append( axes2[0]*np.cos(ang)*np.sin(phi2)+axes2[1]*np.sin(ang)*np.cos(phi2)+center2[1]  )
        zz2.append(raylength(0.128, yy2[-1], rindex_our(0.128)))

    plt.plot(bestx_our, besty_our, "*", ms=8, color='darkviolet', label="best fit(Model 2)")
    #plt.plot(x2, y2, "o", ms=4, color='purple')
    plt.plot(xx2, yy2, "-", color='mediumorchid', label='1sigma allowed region (Model 2)')
    #bestz_our = raylength(0.128, besty_our, rindex_our(0.128))
    #print("Model2 delte = %.3f, rindex@128nm = %.3f, raylength@128nm = %.3f cm" %(besty_our, rindex_our(0.128), raylength(0.128, besty_our, rindex_our(0.128))) )
    #plt.plot(bestx_our, bestz_our, "*", ms=8, color='darkviolet', label="best fit(Model 2)")
    #plt.plot(xx2, zz2, "-", color='mediumorchid', label='1sigma allowed region(Model 2)')
    """

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")

    #plt.xlim(126.45, 126.57)
    #plt.ylim(140.07, 140.17)
    #plt.xlim(50, 70)
    plt.xlim(0.23, 0.30)
    plt.ylim(0.92, 0.96)
    plt.xlabel(r"$\delta$")
    plt.ylabel('Xe absorption ratio')
    #plt.ylabel(r"$L_{Ray}/cm$")
    #plt.xlabel("mu1"); plt.ylabel("mu2");
    plt.legend()
    plt.grid(True)


    #plt.savefig("profile2d_deltaratio_model1.pdf")
    plt.show()
