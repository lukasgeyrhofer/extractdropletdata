#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys,math
from scipy.optimize import curve_fit

import millidrop_dataclass as mdc


def logisticgrowth(t,la0,ln0,lk0,t0):
    return np.exp(ln0 + np.exp(la0)*(t-t0)) / (1.+np.exp(ln0)*(np.exp(np.exp(la0)*(t-t0))-1.)/np.exp(lk0))


def rsquared_log(t,b,parameters):
    ss_res = np.sum((b - logisticgrowth(t,parameters[0],parameters[1],parameters[2],parameters[3]))**2)
    ss_tot = np.sum((b - np.mean(b))**2)
    if ss_tot > 0:  return 1. - ss_res/ss_tot
    else:           return None
    
    
def rsquared_exp(t,b,gr,offset):
    ss_res = np.sum((np.log(b) - gr*t + offset)**2)
    ss_tot = np.sum((np.log(b) - np.mean(np.log(b)))**2)
    if ss_tot > 0:  return 1. - ss_res/ss_tot
    else:           return None
    
    
parser = argparse.ArgumentParser()
parser = mdc.AddCommandLineParameters(parser)

aparser = parser.add_argument_group(description = "==== Algorithm parameters ====")
aparser.add_argument("-m","--maxfev",default=5000,type=int)
aparser.add_argument("-Y","--computeyield",default=False,action="store_true")
aparser.add_argument("-T","--R2threshold",default=None,type=float)
aparser.add_argument("-w","--write_values_to_outfile",default=False,action="store_true")

ffparser = aparser.add_mutually_exclusive_group()
ffparser.add_argument("-E","--exponential", default = False, action = "store_true")
ffparser.add_argument("-L","--logistic",    default = False, action = "store_true")

hparser = parser.add_argument_group(description = "==== Histogram parameters ====")
hparser.add_argument("-b","--bins",default=10,type=int)
hparser.add_argument("-R","--histogramrange",nargs=2,type=float,default=None)

args = parser.parse_args()
data = mdc.DropletData(**vars(args))

if not args.exponential and not args.logistic:
    mode = "exponential"
else:
    mode = "exponential" if args.exponential else "logistic"

growthrates = dict()
yields      = dict()

for experimentLabel, trajectories in data:
    if not growthrates.has_key(experimentLabel):
        growthrates[experimentLabel] = list()
        yields[experimentLabel] = list()
    i = 0
    for trajectory in trajectories:
        t = trajectory[:,0] / args.timerescale
        b = trajectory[:,1]
        
        if mode == "exponential":
            if len(t) >= 2:
                sx  = np.sum(t)
                sxx = np.sum(t*t)
                sy  = np.sum(np.log(b))
                sxy = np.sum(t * np.log(b))
                n   = len(t)
                gr  = (n * sxy - sx * sy)/(n*sxx - sx*sx)
                offset = (gr*sx - sy)/n

                r2 = rsquared_exp(t,b,gr,offset)
                if (args.R2threshold is None) or (args.R2threshold <= r2):
                    growthrates[experimentLabel].append(gr)
                
                if args.verbose:
                    print "{:20s} {:5d} {:12.4e} {:12.4e} {:8.6f}".format(experimentLabel,i,gr,np.exp(offset),r2)
                    i += 1
                
        elif mode == "logistic":
            # length of trajectory has to contain more points than number of fitting parameters
            if len(t) > 4:
                ic = np.array([np.log(abs(np.log(b[0]/b[1])/(t[0] - t[1]))),np.log(b[0]),np.log(b[-1]),t[0]])
                fitMEAN,fitCOV = curve_fit(logisticgrowth,t,b,p0 = ic,maxfev = args.maxfev)
                r2 = rsquared_log(t,b,fitMEAN)
                if (args.R2threshold is None) or (args.R2threshold <= r2):
                    growthrates[experimentLabel].append(np.exp(fitMEAN[0]))
                    yields[experimentLabel].append(np.exp(fitMEAN[2]))
                
                if args.verbose:
                    print "{:20s} {:5d} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:8.6f}".format(experimentLabel,i,np.exp(fitMEAN[0]),np.exp(fitMEAN[1]),np.exp(fitMEAN[2]),fitMEAN[3],r2)
                    i += 1
                

    if not args.verbose:
        print "{:15s} {:.4f} (± {:.4f}) 1/h".format(experimentLabel,np.mean(growthrates[experimentLabel]),np.sqrt(np.std(growthrates[experimentLabel])))

    if len(growthrates[experimentLabel]) > 1:
        if args.histogramrange is None:
            r = (.99 * np.min(growthrates[experimentLabel]), 1.01 * np.max(growthrates[experimentLabel]))
        else:
            r = args.histogramrange
        h,b = np.histogram(growthrates[experimentLabel],bins=args.bins,range = r,density = True)
        b = b[:-1] + np.diff(b)/2.
        outfilename = data.outbasename + experimentLabel + ".growthrates"
        np.savetxt(outfilename,np.transpose([b,h]),fmt = '%.6e')
        
        if args.write_values_to_outfile:
            np.savetxt(data.outbasename + experimentLabel + '.allgrowthrates',np.transpose(growthrates[experimentLabel]),fmt = '%.6e')
    else:
        print "# could not find enough trajectories for histogram"
    if args.computeyield:
        if (len(yields[experimentLabel]) >= 2) and (mode == "logistic"):
            r = (0,1)
            h,b = np.histogram(yields[experimentLabel],bins=args.bins,range = r,density = True)
            b = b[:-1] + np.diff(b)/2.
            outfilename = data.outbasename + experimentLabel + ".yields"
            np.savetxt(outfilename,np.transpose([b,h]),fmt = '%.6e')
            
            if args.write_values_to_outfile:
                np.savetxt(data.outbasename + experimentLabel + '.allyields',yields[experimentLabel], fmt = '%.6e')
            
            



