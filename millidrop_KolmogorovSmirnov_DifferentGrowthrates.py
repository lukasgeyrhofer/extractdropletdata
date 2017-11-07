#!/usr/bin/env python

import numpy as np
import argparse
from scipy.stats import ks_2samp
from termcolor import colored

def col(passed):
    if passed:  return 'red'
    else:       return 'green'

parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*")
parser.add_argument("-t","--pvaluethreshold",default=.05,type=float)
args = parser.parse_args()

data   = list()
labels = list()

for filename in args.infiles:
    try:
        growthrates = np.genfromtxt(filename)
    except:
        continue
    
    data.append(growthrates)
    labels.append(filename.split('.')[0])

nsample = len(data)
D       = np.zeros((nsample,nsample))
pval    = np.zeros((nsample,nsample))

for i in range(nsample):
    for j in range(i):
        D[i,j],pval[i,j] = ks_2samp(data[i],data[j])

print '{:>12s}'.format('') + ' '.join(['{:>12s}'.format(x) for x in labels])
print ''.join(['-' for i in range(13 * len(labels) + 12 )])

for i in range(nsample):
    print '{:>12s}'.format(labels[i]) + ' '.join([colored('{:12.3e}'.format(x),col(c)) for x,c in zip(pval[i],pval[i] < args.pvaluethreshold) if x > 0])
    print '{:>12s}'.format('')        + ' '.join(['{:12.3e}'.format(x) for x,c in zip(D[i],pval[i] < args.pvaluethreshold) if x > 0])
    print




