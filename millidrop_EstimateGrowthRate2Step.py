#!/usr/bin/env python

import numpy as np
import argparse
import pandas as pd
import sys,math
from scipy.optimize import curve_fit

import millidrop_dataclass_new as mdc


# sigmoid with minimum value.
# a few parameters should be postive, thus takes their square here
# furthermore, fit on logarithm of sigmoid curve, this should put more emphasis on initial times
def sigmoid(t,rate,centerpos,minval2,maxval2):
    return np.log(minval2**2 + (maxval2**2 - minval2**2)/(1 + np.exp(- rate*(t-centerpos))))

def extend(x,fulllength):
    start = int((len(x)-fulllength)/2)
    return np.concatenate([np.zeros(start),x,np.zeros(fulllength-len(x)-start)])

def compute_weights(growthrates,numberpoints):
    w = np.exp(-np.concatenate([np.ones(numberpoints/2)*growthrates[0],growthrates,np.ones(numberpoints/2)*growthrates[-1]]))
    w /= np.sum(w)
    return w

parser = argparse.ArgumentParser()
parser_data = parser.add_argument_group(description = "==== Data parameters ====")
parser_data.add_argument("-i", "--infiles"          )
parser_data.add_argument("-t", "--templatefile",    default=None)
parser_data.add_argument("-u", "--timerescale", default=3.6e3, type=float)
parser_data.add_argument("-E", "--excludelabels",default=None,nargs ="*")

parser_gr = parser.add_argument_group(description = "==== Growthrate Estimation ====")
parser_gr.add_argument("-n","--numberpoints",default=4,type=int)
parser_gr.add_argument("-f","--upperfrac",default=.95,type=float)
parser_gr.add_argument("-M","--maxfev",default=1000,type=int)

parser_io = parser.add_argument_group(description = "==== I/O ====")
parser_io.add_argument("-o","--outbasename",     default=None)
parser_io.add_argument("-S","--skipsecond",default=False,action="store_true")
parser_io.add_argument("-H","--compute_histogram",default=False,action="store_true")
parser_io.add_argument("-R","--histo_range",default=[0,2],nargs=2,type=float)
parser_io.add_argument("-B","--histo_bins",default=40,type=int)
args = parser.parse_args()

data = mdc.DropletDataNew(**vars(args))
if len(args.excludelabels) > 0:
    for el in args.excludelabels:
        data.ExcludeLabel(el)

for label in data.labels:
    gr_all = list()
    gr_all_2step = list()
    for dropdata in data[label]:
        # dropdata contains the whole pandas object loaded from the csv-file
        # separate into back/forth, choose channel of signal, rescale time
        # get two np.arrays back
        t1,t2 = mdc.getTrajectories(dropdata,'fluo_3_median',timerescale = args.timerescale)
        
        # first step:
        # compute growthrates as 'local' average of exponential growthrates over 'args.numberpoints' points
        gr1 = np.array([mdc.MLSQ_fit(t1[i:i+args.numberpoints,0],np.log(t1[i:i+args.numberpoints,1]),returnoffset = False) for i in range(len(t1[:,0])-args.numberpoints)],dtype=np.float)
        gr2 = np.array([mdc.MLSQ_fit(t2[i:i+args.numberpoints,0],np.log(t2[i:i+args.numberpoints,1]),returnoffset = False) for i in range(len(t2[:,0])-args.numberpoints)],dtype=np.float)

        #np.savetxt(args.outbasename + '{:04d}-{}F'.format(i,label),np.transpose([t1[:,0],t1[:,1]]))
        #np.savetxt(args.outbasename + '{:04d}-{}B'.format(i,label),np.transpose([t2[:,0],t2[:,1]]))
        

        # add to list
        gr_all.append(mdc.MeanUpperFrac(gr1,args.upperfrac))
        gr_all.append(mdc.MeanUpperFrac(gr2,args.upperfrac))

        if not args.skipsecond:
            # second step:
            # use growthrates as inverse weights to fit global sigmoid
            
            # estimate start values for fit
            # a few values need to be restricted to positive values, thus they are squared in the function above, need to take root here
            param0 = np.array([(t1[:,0])[((t1[:,1]-t1[0,1])*(t1[-1,1]-t1[:,1])).argmax()],1,np.sqrt(t1[0,1]),np.sqrt(t1[-1,1])])
            
            # define weights as inverse growthrates
            weights1 = compute_weights(gr1,args.numberpoints)
            weights2 = compute_weights(gr2,args.numberpoints)
            
            # fit sigmoids
            try:
                fit1,cov1 = curve_fit(sigmoid,t1[:,0],np.log(t1[:,1]),p0 = param0,sigma = weights1,maxfev = args.maxfev)
                gr_all_2step.append(fit1[0])
            except:
                pass
            
            try:
                fit1,cov1 = curve_fit(sigmoid,t1[:,0],np.log(t1[:,1]),p0 = param0,sigma = weights1,maxfev = args.maxfev)
                gr_all_2step.append(fit1[0])
            except:
                pass
            
    np.savetxt(args.outbasename + label + '_growthrates',gr_all)
    if args.compute_histogram:
        h,b = np.histogram(gr_all,range = args.histo_range,bins = args.histo_bins, density = True)
        b = .5*(b[1:] + b[:-1])
        np.savetxt(args.outbasename + label + '_histogram',np.transpose([b,h]))
    
    if not args.skipsecond:
        np.savetxt(args.outbasename + label + '_growthrates_2step',gr_all_2step)
        if args.compute_histogram:
            h,b = np.histogram(gr_all_2step,range = args.histo_range,bins = args.histo_bins, density = True)
            b = .5*(b[1:] + b[:-1])
            np.savetxt(args.outbasename + label + '_histogram_2step',np.transpose([b,h]))
