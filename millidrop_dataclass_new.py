#!/usr/bin/env python

import numpy as np
import argparse
import pandas as pd
import sys,math
import os


class DropletDataNew(object):
    def __init__(self,**kwargs):
        self.__dropdir                    = kwargs.get("infiles",None)
        # templatefile to assign the labels from different wells
        self.__templatefile               = kwargs.get("templatefile",None)

        # downstream processing, not actively accessed by this data object
        # used in several of the derivative scripts that use this object, only used to store this value
        self.__timerescale                = kwargs.get("timerescale",3.6e3)
        # store folder to output files
        self.__outbasename                = kwargs.get("outbasename","")
        if self.__outbasename is None:
            self.__outbasename = ""
        
        # further options when loading the data
        self.__splitBackForthTrajectories = kwargs.get("SplitBackForthTrajectories",True)
        
        self.__excludelabels = list()
        
        self.__dropmap  = self.LoadDropMap(self.__templatefile)
        self.__idropmap = self.InvertDropmap(self.__dropmap)
        
        self.__data     = dict()
        for i in range(len(self.__dropmap)):
            fn = os.path.join(self.__dropdir,'{:04d}.csv'.format(i))
            self.__data[i] = pd.read_csv(fn)
            if 'time' in self.__data[i]:
                self.__data[i]['time'] /= self.__timerescale

    def LoadDropMap(self,filename):
        tpfile = pd.read_csv(filename,header = 4)
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

    def InvertDropmap(self,dropmap):
        idm = dict()
        for i,wl in enumerate(dropmap):
            if not idm.has_key(wl[1]):
                idm[wl[1]] = list()
            idm[wl[1]].append(i)
        for keys in idm.iterkeys():
            idm[keys] = np.array(idm[keys],dtype=int)
        return idm
    
    def ExcludeLabel(self,label):
        self.__excludelabels.append(label)
    
    def __getattr__(self,key):
        if key == 'labels':
            return [label for label in self.__idropmap.keys() if label not in self.__excludelabels]
        #else:
            #super(DropletDataNew,self).__getattr__(key)
            
    
    def __getitem__(self,key):
        if key in self.__idropmap.keys():
            for dropID in self.__idropmap[key]:
                yield self.__data[dropID]
        else:
            raise KeyError('"{}" not in allowed labels: ('.format(key) + ', '.join(self.__idropmap.keys()) + ')')


def getTrajectories(data,channel,SplitBackForthTrajectories=True,timerescale = 1):
    if channel in data.columns and 'time' in data.columns:
        t = np.array(data['time']/timerescale,dtype=np.float)
        x = np.array(data[channel],dtype=np.float)
        
        a = np.transpose([t,x])
        
        if SplitBackForthTrajectories:
            return a[::2,:],a[1::2,:]
        else:
            return a


    else:
        return None


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

def MeanUpperFrac(x,upperfrac = .95,returnnumber = False):
    if len(x) > 0:
        uf_val = np.max(x) * upperfrac
        m = np.mean(x[x >= uf_val])
        n = len(x[x >= uf_val])
        if returnnumber:
            return m,n
        else:
            return m
    else:
        return None

# helper script to add all cmdline parameters
def AddCommandLineParameters(parser):
    ioparser = parser.add_argument_group(description = "==== I/O parameters ====")
    
    ioparser.add_argument("-i", "--infiles"          )
    ioparser.add_argument("-t", "--templatefile",    default=None)
    ioparser.add_argument("-r", "--restrictionfile", default=None)
    ioparser.add_argument("-o", "--outbasename",     default=None)
    
    ioparser.add_argument("-C", "--datacolumns", nargs="*",type=str)
    ioparser.add_argument("-u", "--timerescale", default=3.6e3, type=float)
    
    ioparser.add_argument("-B", "--SplitBackForthTrajectories", default = False, action = "store_true")
    ioparser.add_argument("-H", "--NonHiccupLoading",           default = False, action = "store_true")
    ioparser.add_argument("-D", "--IgnoreAdditionalDroplets",   default = False, action = "store_true")
    
    ioparser.add_argument("-v","--verbose",                     default = False, action = "store_true")
    
    return parser


