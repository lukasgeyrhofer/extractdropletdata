#!/usr/bin/env python

import argparse
import numpy as np
import sys,math
from scipy.stats import poisson


def growthrate(time,n0,n1):
    return np.log(n1/n0)/time

def diff(x):
    return 0.5*np.diff(np.concatenate([np.array([0]),x]) + np.concatenate([x,np.array([1])]))


parser = argparse.ArgumentParser()
parser.add_argument("-n","--initialconditionsPoisson",default=25,type=int)
parser.add_argument("-N","--Nmax",default=100,type=int)
parser.add_argument("-t","--thresholdtimefile",default=None)
parser.add_argument("-T","--threshold",type=float,default=.5)
parser.add_argument("-S","--signaltocellnumber",type=float,default=1e4) # this might not be completely correct

parser.add_argument("-a","--amin",default=0,type=float)
parser.add_argument("-A","--amax",default=1,type=float)
parser.add_argument("-d","--da",default=.01,type=float)
args = parser.parse_args()

try:
    thdata = np.genfromtxt(args.thresholdtimefile)
except:
    raise IOError

tlist  = thdata[:,0]
dt     = tlist[1] - tlist[0]
ptlist = thdata[:,1] * dt

if args.Nmax > 0:
    nlist  = np.arange(1,args.Nmax)
    pnlist = poisson.pmf(nlist,args.initialconditionsPoisson)
else:
    nlist  = np.array([args.initialconditionsPoisson])
    pnlist = np.ones(1)

alist  = np.arange(start = args.amin,stop = args.amax, step = args.da)
Calist = np.zeros(len(alist))

for i in range(len(alist)):
    for n,pn in zip(nlist,pnlist):
        for t,pt in zip(tlist,ptlist):
            if growthrate(t,n,args.threshold*args.signaltocellnumber) <= alist[i]:
                Calist[i] += pn*pt

Palist  = diff(Calist)
Palist /= np.sum(Palist)*args.da

for a,pa in zip(alist,Palist):
    print "{:.3f} {:.5e}".format(a,pa)

