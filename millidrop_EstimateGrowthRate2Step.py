#!/usr/bin/env python

import numpy as np
import argparse
import pandas as pd
import sys,math
from scipy.optimize import curve_fit

import millidrop_dataclass_new as mdc


# sigmoid with minimum value.
# a few parameters should be postive, thus takes their square here
def sigmoid(t,center,width2,minval2,maxval2):
    return minval2**2 + (maxval2**2 - minval2**2)/(1 + np.exp(-(t-center)/width2**2))

parser = argparse.ArgumentParser()
parser = mdc.AddCommandLineParameters(parser)
parser_gr = parser.add_argument_group(description = "==== Growthrate Estimation ====")
parser_gr.add_argument("-n","--numberpoints",default=4,type=int)
parser_gr.add_argument("-f","--upperfrac",default=.95,type=float)
parser_gr.add_argument("-M","--maxfev",default=1000,type=int)
args = parser.parse_args()

data = mdc.DropletDataNew(**vars(args))
data.ExcludeLabel('Empty')


for label in data.labels:
    grall = list()
    grall2step = list()
    for dropdata in data[label]:
        # dropdata contains the whole pandas object loaded from the csv-file
        # separate into back/forth, choose channel of signal, rescale time
        # get two np.arrays back
        t1,t2 = mdc.getTrajectories(dropdata,'fluo_3_median',timerescale = args.timerescale)
        
        # first step:
        # compute growthrates as 'local' average of exponential growthrates over 'args.numberpoints' points
        gr1 = np.array([mdc.MLSQ_fit(t1[i:i+args.numberpoints,0],np.log(t1[i:i+args.numberpoints,1]),returnoffset = False) for i in range(len(t1[:,0])-args.numberpoints)],dtype=np.float)
        gr2 = np.array([mdc.MLSQ_fit(t2[i:i+args.numberpoints,0],np.log(t2[i:i+args.numberpoints,1]),returnoffset = False) for i in range(len(t2[:,0])-args.numberpoints)],dtype=np.float)

        # add to list
        grall.append(mdc.MeanUpperFrac(gr1,args.upperfrac))
        grall.append(mdc.MeanUpperFrac(gr2,args.upperfrac))


        # second step:
        # use growthrates as inverse weights to fit global sigmoid
        
        # estimate start values for fit
        # a few values need to be restricted to positive values, thus they are squared in the function above, need to take root here
        param0 = np.array([(t1[:,0])[((t1[:,1]-t1[0,1])*(t1[-1,1]-t1[:,1])).argmax()],1,np.sqrt(t1[0,1]),np.sqrt(t1[-1,1])])
        
        # define weights as inverse growthrates
        # rescale properly first, and add a small number to avoid dividing by 0
        weights1 = 1./(gr1 - np.min(gr1) + 1e-10)
        weights2 = 1./(gr2 - np.min(gr2) + 1e-10)
        startindex = int(args.numberpoints/2)
        
        # fit sigmoids
        try:
            fit1,cov1 = curve_fit(sigmoid,t1[startindex:len(t1[:,0])-startindex,0],t1[startindex:len(t1[:,0])-startindex,1],p0 = param0,sigma = weights1,maxfev = args.maxfev)
            grall2step.append(1./fit1[0]**2)
        except:
            continue
        
        try:
            fit2,cov2 = curve_fit(sigmoid,t2[startindex:len(t2[:,0])-startindex,0],t2[startindex:len(t2[:,0])-startindex,1],p0 = param0,sigma = weights2,maxfev = args.maxfev)
            grall2step.append(1./fit2[0]**2)
        except:
            continue
        
        # add to list
        
        #fp = open(label + '_growthrates','w')
    np.savetxt(label + '_growthrates',grall)
    np.savetxt(label + '_growthrates_2step',grall2step)
    #print grall
    #print grall2step
