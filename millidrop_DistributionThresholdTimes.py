#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import argparse
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
parser = mdc.AddCommandLineParameters(parser)

aparser = parser.add_argument_group(description = "==== algorithm parameters ====")
aparser.add_argument("-T","--threshold",type=float,default=.4)

hparser = parser.add_argument_group(description = "==== Histogram parameters ====")
hparser.add_argument("-b","--bins",default=10,type=int)
hparser.add_argument("-R","--histogramrange",nargs=2,type=float,default=None)


args = parser.parse_args()
data = mdc.DropletData(**vars(args))

for ExperimentLabel,Trajectories in data:
    thresholdtimes = list()
    for traj in Trajectories:
        t = traj[:,0] / args.timerescale
        b = traj[:,1]
        
        timesabovethreshold = t[b>=args.threshold]
        if len(timesabovethreshold) > 0:
            thresholdtimes.append(np.min(timesabovethreshold))
        
    if len(thresholdtimes) > 0:
        hist,bins = np.histogram(thresholdtimes,bins = args.bins,range = args.histogramrange,density = True)
        bins = bins[:-1] + np.diff(bins)/2.
        
        if args.verbose:
            print "{:15s}: storing histogram of {:d} values, ({:.2f} ± {:.2f})".format(ExperimentLabel,len(thresholdtimes),np.mean(thresholdtimes),np.std(thresholdtimes))
        
        outfilename = data.outbasename + ExperimentLabel + ".threshold{:0.3e}".format(args.threshold)
        np.savetxt(outfilename,np.transpose(np.array([bins,hist])),fmt = "%.6e")


                    







