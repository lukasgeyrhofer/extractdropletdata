#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse
import sys,math
import os

import pylab as p

import millidrop_dataclass as mdc

def LoadDropMap(filename):
    try:
        tpfile = pd.read_csv(filename,header = 4)
    except:
        raise IOError("could not open templatefile")

    tpFileOrd=tpfile.set_index('order')
    dropMap = list()
    for i in range(0,len(tpFileOrd)-1,2):
        for j in range(2*(tpFileOrd.droplet_number[i])):
            if j%2==0 :
                dropMap.append([tpFileOrd.well[i], tpFileOrd.description[i]])
            else :
                dropMap.append([tpFileOrd.well[i+1], tpFileOrd.description[i+1]])
    return np.array(dropMap)

def InvertDropmap(dropmap):
    idm = dict()
    for i,wl in enumerate(dropmap):
        if not idm.has_key(wl[1]):
            idm[wl[1]] = list()
        idm[wl[1]].append(i)
    for keys in idm.iterkeys():
        idm[keys] = np.array(idm[keys],dtype=int)
    return idm

def MLSQ_fit(x,y,returnoffset = True):
    assert len(x) == len(y)
    sx = np.sum(x)
    sxx = np.dot(x,x)
    sy = np.sum(y)
    sxy = np.dot(x,y)
    n = len(x)
    a = (n*sxy - sx*sy)/(n*sxx-sx*sx)
    b = (sy-a*sx)/n
    if returnoffset:
        return a,b
    else:
        return a

def FracMean(x,frac = .95,returnnumber = False):
    mx = np.max(x)
    m  = np.mean(x[x >= frac * mx])
    if returnnumber:
        n = len(x[x >= frax * mx])
        return m,n
    else:
        return m

def Statistics(x):
    return np.mean(x),np.std(x),np.median(x)

def hist_s(x,range = (0,1),bins =20):
    h,b = np.histogram(x,range=range,bins=bins)
    b = b[:-1] + np.diff(b)
    return h,b

parser = argparse.ArgumentParser()
parser.add_argument("-d","--dropletDir",required=True)
parser.add_argument("-t","--template",required=True)

parser.add_argument("-n","--numberpoints",default=5,type=int)
parser.add_argument("-f","--maxfrac",default=.95,type=float)
args = parser.parse_args()

dropmap  = LoadDropMap(args.template)
idropmap = InvertDropmap(dropmap)
del idropmap['Empty']

df = dict()
for i,fn in enumerate([args.dropletDir + f for f in os.listdir(args.dropletDir) if os.path.isfile(os.path.join(args.dropletDir,f))]):
    if i<=len(dropmap):
        df[i] = pd.read_csv(fn)

for label in idropmap.iterkeys():
    grlist = list()
    for i in idropmap[label]:
        t = df[i]['time'] / 3.6e3
        f = df[i]['fluo_1_mean']
        if len(t) >= args.numberpoints:
            gr = np.array( [MLSQ_fit(t[i:i+args.numberpoints],np.log(f[i:i+args.numberpoints]),returnoffset=False) for i in range(len(t)-args.numberpoints)], dtype=np.float)
            grlist.append(FracMean(gr,frac = args.maxfrac))
    #print "{:20s} {:.6f} {:.6f} {:.6f}".format(label,*Statistics(grlist))
    
    print grlist
    
        
