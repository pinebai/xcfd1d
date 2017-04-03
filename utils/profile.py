import sys
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) == 1:
    varnum = 1
else:
    varnum = int(sys.argv[1])

line = []
varname = []
varsymb = []
varunit = []

with open('output.dat') as fh:
    nvar = int(fh.readline())
    for i in range(nvar):
        buffer = fh.readline()
        line = buffer.split(',')
        varname.append(line[0])
        varsymb.append(line[1])
        varunit.append(line[2])
    nexa = int(fh.readline())
    dexa = np.zeros( (nvar,nexa) )
    for i in range(nexa):
        buffer = fh.readline()
        line = buffer.split()
        for j in range(nvar):
            dexa[j,i] = float(line[j])
    ndat = int(fh.readline())
    data = np.zeros( (nvar,ndat) )
    for i in range(ndat):
        buffer = fh.readline()
        line = buffer.split()
        for j in range(nvar):
            data[j,i] = float(line[j])

plt.plot(dexa[0,:], dexa[varnum,:], color='k', label='Exact')
plt.plot(data[0,:], data[varnum,:], color='g', label='Computed')

xlabel = varname[0] + ' ' + varunit[0]
ylabel = varname[varnum] + ' ' + varunit[varnum]

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.grid(True)
plt.legend(loc='upper right')
plt.show()
