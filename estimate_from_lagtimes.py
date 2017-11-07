#!/usr/bin/env python
# -*- coding: utf-8 -*-

# estimate value for tau from Kjeldsen et al 2014

import argparse
import numpy as np
from scipy.optimize import curve_fit
import sys,math


def lagtime(conc,tau,mic):
    def k(conc1,mic):
        return params['lambda'] * np.log(conc1/mic)
    return (1 + 1./k(conc,mic)) * (1+np.power(k(conc,mic),3)/(tau*tau*params['lambda']*params['inoculum']))


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile")
parser.add_argument("-l","--ABlambda",default=5,type=float)
parser.add_argument("-n","--inoculum",default=1e5,type=float)
parser.add_argument("-M","--maxfev",default=5000,type=int)
args = parser.parse_args()

global params
params = {'lambda':args.ABlambda,'inoculum':args.inoculum}


try:
    data = np.genfromtxt(args.infile)
except:
    raise IOError


p0 = [1.,1.]
np.seterr(all = 'ignore')
a,b = curve_fit(lagtime,data[1:,0],data[1:,1]-data[0,1],p0=p0,maxfev = args.maxfev)

print "tau    {:.3e} (± {:.3e})".format(abs(a[0]),np.sqrt(b[0,0]))
print "mic    {:.3e} (± {:.3e})".format(abs(a[1]),np.sqrt(b[1,1]))
