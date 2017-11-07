#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
parser = mdc.AddCommandLineParameters(parser)


parser.add_argument("-S","--lowerSignalCutoff",default=None,type=float)
parser.add_argument("-R","--histogramRange",default = None,nargs=2,type=float)
parser.add_argument("-b","--histogramBins",default=10,type=int)
args = parser.parse_args()

data = mdc.DropletData(**vars(args))

hRange = args.histogramRange
hBins  = args.histogramBins


for experimentlabel,trajectories in data:
    finalSignal = list()
    for trajectory in trajectories:
        if not args.lowerSignalCutoff is None:
            if len(trajectory) > 1:
                if trajectory[-1][1] > args.lowerSignalCutoff:
                    finalSignal.append(trajectory[-1][1])
        else:
            if len(trajectory) > 1:
                finalSignal.append(trajectory[-1][1])
    
    if len(finalSignal) > 1:
        outfilename = data.outbasename + experimentlabel + ".yield"
        h,b = np.histogram(finalSignal,bins = hBins,range = hRange)
        b = b[:-1] + .5 * np.diff(b)
        np.savetxt(outfilename,np.transpose([b,h]))
        
        stat = [np.mean(finalSignal),np.std(finalSignal)]
    else:
        stat = [np.nan,np.nan]
    
    if args.verbose:
        print "{:15s} {:.4f} (± {:.4f}) from {:d} droplets".format(experimentlabel,stat[0],stat[1],len(finalSignal))





