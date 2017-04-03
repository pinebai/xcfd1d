import sys
import matplotlib.pyplot as plt
import numpy as np

nvar = 11
ndat = sum(1 for line in open('residual.dat'))
ndat -= 1

line = []

with open('residual.dat') as fh:
    buffer = fh.readline()
    data = np.zeros( (nvar,ndat) )
    for i in range(ndat):
        buffer = fh.readline()
        line = buffer.split()
        for j in range(nvar):
            data[j,i] = float(line[j])

plt.plot(data[0,:], data[2,:], color='r', label='rho')
plt.plot(data[0,:], data[3,:], color='g', label='du')
plt.plot(data[0,:], data[4,:], color='b', label='E')

xlabel = 'Number of Iterations'
ylabel = 'L1 Residual'

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()
