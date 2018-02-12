import numpy as np
import math
import time
import re
import sys, getopt
import matplotlib.pyplot as plt

infile=str(sys.argv[1])
outfile=str(sys.argv[2])
deleps=float(sys.argv[3])


mineps=-10.0
maxeps=10.0

with open(infile, 'r') as f:
    fft = f.read()  # read the whole file
ffl= fft.splitlines()  # split by lines
nlines = len(ffl)  # number of lines
fea=[]
for line in ffl :
    sline=line.split()

    fea.append(float(sline[1]))

fea=np.array(fea)

epsl=np.linspace(mineps,maxeps,np.floor((maxeps-mineps)/deleps+1))

saveff=[]
for eps in epsl:
    ff=-np.log(np.sum(np.exp(-fea-eps*(range(nlines)))))
    saveff.append((eps,ff))




with open(outfile, 'w') as f:
    for ii in saveff:
        f.write(str(ii[0])+"  " +str(ii[1]) + "\n")

#plt.plot(saveff)
#plt.show()
