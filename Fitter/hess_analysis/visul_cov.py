#!/usr/bin/env python
# coding=utf-8

import numpy as np

def load_cov(filename):
    with open(filename) as f:
        cov = []
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            tmp = []
            for elem in data:
                tmp.append(float(elem))
            cov.append(tmp)

    cov = np.array(cov)
    return cov


import matplotlib.pyplot as plt
import matplotlib.colors as colors

def main():
    cov = load_cov("hessian_matrix.txt")

    XX = ["a0", "aUV", "aIR", "delta", "ratio", "mu1", "sigma1", "A2", "mu2", "sigma2", "nu_f"]
    YY = ["a0", "aUV", "aIR", "delta", "ratio", "mu1", "sigma1", "A2", "mu2", "sigma2", "nu_f"]

    #cov_part = cov[:4, :4]
    cov_part = cov
    fig, ax = plt.subplots()
    im = ax.imshow(cov_part)  
        #norm=colors.LogNorm() )

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(XX)))
    ax.set_yticks(np.arange(len(YY)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(XX)
    ax.set_yticklabels(YY) 

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

    cbar = ax.figure.colorbar(im, ax=ax)

    #for i in range(len(XX)):
    #    for j in range(len(YY)):
    #        text = ax.text(j, i, cov_part[i, j],
    #            ha="center", va="center", color="w")

    ax.set_title("Hessian Matrix")
    fig.tight_layout()
    plt.savefig("hesMatrix.pdf")

    plt.show()

if __name__ == '__main__':
    main()

