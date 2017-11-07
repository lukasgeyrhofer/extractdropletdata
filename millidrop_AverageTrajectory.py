#!/usr/bin/env python

import argparse
import numpy as np
import sys,math

import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
parser = mdc.AddCommandLineParameters(parser)

tparser = parser.add_argument_group(description = "==== New grid parameters ====")
tparser.add_argument("-M","--maxtime",default=None,type=float)
tparser.add_argument("-m","--mintime",default=0,type=float)
tparser.add_argument("-n","--datapoints",default=50,type=int)

aparser = parser.add_argument_group(description = "==== Algorithm parameters ====")
aparser.add_argument("-d","--stddev",default=False,action="store_true")

args = parser.parse_args()
data = mdc.DropletData(**vars(args))


for label,trajectories in data:
    n = len(trajectories)
    if n >= 1:
        if args.maxtime is None:
            maxtime = np.max([t[-1,0] for t in trajectories if len(t) >= 1])
        else:
            maxtime = args.maxtime
        timegrid        = np.linspace(start = args.mintime,stop = maxtime,num = args.datapoints)
        sumtrajectories = np.zeros((np.shape(trajectories[0])[1]-1,args.datapoints))
        if args.stddev:
            sum2trajectories = np.zeros((np.shape(trajectories[0])[1]-1,args.datapoints))
        for t in trajectories:
            if len(t) > 0:
                for i in range(1,np.shape(trajectories[0])[1]):
                    values = np.interp(timegrid,t[:,0],t[:,i])
                    sumtrajectories[i-1] += values
                    if args.stddev:
                        sum2trajectories[i-1] += values*values
            else:
                n -= 1
    if n >= 1:
        outfilename = data.outbasename + label + ".average"
        print "{:12s}: saving average from {:d} trajectories to file '{:s}'".format(label,n,outfilename)
        outdata = np.array([timegrid/args.timerescale])
        outdata = np.concatenate([outdata,np.array(sumtrajectories/n)],axis=0)
        np.savetxt(outfilename,np.transpose(outdata),fmt = '%.6e')
        if args.stddev:
            outfilename = outbasename + label + ".stddev"
            outdatastddev = np.array([timegrid/args.timerescale])
            outdatastddev = np.concatenate([outdatastddev,np.sqrt(n*sum2trajectories - sumtrajectories*sumtrajectories)/np.sqrt(n*n-n)],axis=0)
            np.savetxt(outfilename,np.transpose(outdatastddev),fmt = '%.6e')
            
        

