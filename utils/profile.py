import sys
import matplotlib.pyplot as plt
import numpy as np

colrs = ['b', 'g', 'r', 'c', 'm', 'y']

filenm = []

if len(sys.argv) >= 6:
    print 'Maximum ten profiles allowed'
    sys.exit()
elif len(sys.argv) >= 3:
    filnum = len(sys.argv)-2
    varnum = int(sys.argv[1])
    for j in range(2,len(sys.argv)):
        filenm.append(sys.argv[j])
else:
    print 'Arguments not correct'
    sys.exit()

line = []
varname = []
varsymb = []
varunit = []

for j in range(filnum):
    with open(filenm[j]) as fh:
        nvar = int(fh.readline())
        for i in range(nvar):
            buffer = fh.readline()
            line = buffer.split(',')
            varname.append(line[0])
            varsymb.append(line[1])
            varunit.append(line[2])

        nexa = int(fh.readline())
        dexa = np.zeros((2,nexa))
        for i in range(nexa):
            buffer = fh.readline()
            line = buffer.split()
            dexa[0,i] = float(line[0])
            dexa[1,i] = float(line[varnum])
        if j == 0:
            plt.plot(dexa[0,:], dexa[1,:], color='k', label='Exact')

        ndat = int(fh.readline())
        data = np.zeros((2,ndat))
        for i in range(ndat):
            buffer = fh.readline()
            line = buffer.split()
            data[0,i] = float(line[0])
            data[1,i] = float(line[varnum])

        plt.plot(data[0,:], data[1,:], color=colrs[j], label='Computed')

xlabel = varname[0] + ' ' + varunit[0]
ylabel = varname[varnum] + ' ' + varunit[varnum]

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.grid(True)
plt.legend(loc='upper right')
plt.show()
