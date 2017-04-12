import sys
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) == 3:
    var = sys.argv[1]
    fnm = sys.argv[2]
else:
    print 'Arguments not correct'
    sys.exit()

nvar = 11
ndat = sum(1 for line in open(fnm))
ndat -= 1

line = []

with open(fnm) as fh:
    buffer = fh.readline()
    data = np.zeros( (nvar,ndat) )
    for i in range(ndat):
        buffer = fh.readline()
        line = buffer.split()
        for j in range(nvar):
            data[j,i] = float(line[j])

xlabel = 'Number of Iterations'

if var == 'l1':
    plt.plot(data[0,:], data[2,:], color='r', label='rho')
    plt.plot(data[0,:], data[3,:], color='g', label='du')
    plt.plot(data[0,:], data[4,:], color='b', label='E')
    ylabel = 'L1 Residual'
elif var == 'l2':
    plt.plot(data[0,:], data[5,:], color='r', label='rho')
    plt.plot(data[0,:], data[6,:], color='g', label='du')
    plt.plot(data[0,:], data[7,:], color='b', label='E')
    ylabel = 'L2 Residual'
else:
    plt.plot(data[0,:], data[8,:], color='r', label='rho')
    plt.plot(data[0,:], data[9,:], color='g', label='du')
    plt.plot(data[0,:], data[10,:], color='b', label='E')
    ylabel = 'Max Residual'

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()
