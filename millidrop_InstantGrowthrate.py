#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
parser = mdc.AddCommandLineParameters(parser)

aparser = parser.add_argument_group(description = "==== Algorithm parameters ====")

hparser = parser.add_argument_group(description = "==== Histogram parameters ====")
hparser.add_argument("-B","--bins",default=120,type=int)
hparser.add_argument("-R","--histogramrange",nargs=2,type=float,default=[-.1,1.1])

args = parser.parse_args()
data = mdc.DropletData(**vars(args))


for experimentLabel, trajectories in data:
    instantgrowthrate = list()
    for trajectory in trajectories:
        for i in range(1,len(trajectory)):
            instantgrowthrate.append(np.log(trajectory[i,1]/trajectory[i-1,1])/(trajectory[i,0]-trajectory[i-1,0])*1e3)
    instantgrowthrate = np.array(instantgrowthrate)
    if args.histogramrange is None:
        r = (.99 * np.min(instantgrowthrate), 1.01 * np.max(instantgrowthrate))
    else:
        r = args.histogramrange
    h,b = np.histogram(instantgrowthrate,bins=args.bins,range = r,density = True)
    b = b[:-1] + np.diff(b)/2.
    outfilename = data.outbasename + experimentLabel + ".instantgrowthrates"
    np.savetxt(outfilename,np.transpose([b,h]),fmt = '%.6e')
    

