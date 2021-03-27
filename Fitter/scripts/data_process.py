import numpy as np
import matplotlib.pyplot as plt

index = 0
wl, mean, err = [], [], []
tmp_mean, tmp_up, tmp_low = 0, 0, 0
wl0, wl1, wl2 = 0, 0, 0
with open("../data/tmp.txt") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        if index == 3:
            index = 0
        if index == 0:
            wl0 = float(data[0])
            tmp_mean = float(data[1])
        elif index == 1:
            wl1 = float(data[0])    
            tmp_up = float(data[1])
        elif index == 2:
            wl2 = float(data[0])    
            tmp_low = float(data[1])

        index+=1

        if index == 3:
            wl.append((wl0+wl1+wl2)/3)
            mean.append((tmp_mean + tmp_up + tmp_low)/3.)
            err.append(max(abs(tmp_up-tmp_mean), abs(tmp_mean-tmp_low)))


wl = np.array(wl)
print(wl)
mean = np.array(mean)
err = np.array(err)

with open("../data/cell116mm.txt", "w") as f:
    for i, j, k in zip(wl, mean, err):
        f.write("%.3f %.3f %.3f" %(i, j, k) )
        f.write("\n")

#plt.errorbar(wl, mean, yerr=err)
#plt.show()
