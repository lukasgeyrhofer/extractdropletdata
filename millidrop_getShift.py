#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infiles",nargs="*") # can be multiple files, or also wildcards ...
parser.add_argument("-T","--lowerthreshold",type=float,default=3e-2)
args = parser.parse_args()


for filename in args.infiles:
    try:
        # load data to numpy array
        data = np.genfromtxt(filename,skip_header = 1,delimiter = ',', dtype = float)

        dropletnr = int(filename.replace(".csv",""))
        
        # get data for time and bacterial density
        t = data[:,3]
        b = data[:,11]
        
        # restrict data to remove data below some observation threshold
        t = t[b>args.lowerthreshold]
        b = b[b>args.lowerthreshold]
                
        # generate trajectories from alternating points
        t1 = t[0::2]
        t2 = t[1::2]
        b1 = b[0::2]
        b2 = b[1::2]
        
        # interpolate time to estimate the timestamp if b2 would be on the (t1,b1) trajectory
        tinterp =  np.interp(b2,b1,t1)
        
        # average over all such obtained shifts
        averageshift = np.mean(tinterp - t2)
        
        print dropletnr,averageshift
    except:
        # just skip to next file if something does not work
        continue
    
