#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
args = parser.parse_args()

data = np.array(())
first = True
for filename in args.infiles:
    if first:
        data = np.array([np.genfromtxt(filename)])
        first = False
    else:
        data = np.concatenate((data,[np.genfromtxt(filename)]))

x = np.arange(data.shape[0])

for i in range(data.shape[1]):
    if i%2 == 0:
        time = data[:,i,0]
        # reverse timestamp
        #data[:,i,0] = time[::-1]
        # add fit
        data[:,i,0] = time - (x-500.) * 3.7

for i in range(len(args.infiles)):
    np.savetxt(args.infiles[i].replace("data","revT"),data[i],fmt = '%.6e')




