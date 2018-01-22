#!/usr/bin/env python

import numpy as np
import argparse
import pandas as pd
import sys,math
import os



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
        
        self.__maxtime        = kwargs.get("maxtime",None)
        self.__fitpoints      = kwargs.get("fitpoints",5)
        self.__timecolumnname = kwargs.get("timecolumnname",'time')

        # further options when loading the data
        self.__splitBackForthTrajectories = kwargs.get("SplitBackForthTrajectories",True)
        
        if kwargs.has_key("excludelabels"):
            self.__excludelabels = [l.upper() for l in kwargs.get("excludelabels")]
        else:
            self.__excludelabels = list()
        
        self.__dropmap  = self.LoadDropMap(self.__templatefile)
        self.__idropmap = self.InvertDropmap(self.__dropmap)
        
        self.__columns = None
        
        
        self.__data     = dict()
        for i in range(len(self.__dropmap)):
            fn = os.path.join(self.__dropdir,'{:04d}.csv'.format(i))
            self.__data[i] = pd.read_csv(fn)
            if self.__columns is None:
                self.__columns = list(self.__data[i].columns)
            else:
                for c in self.__data[i].columns:
                    if not c in self.__columns:
                        self.__columns.drop(c)
            if self.__timecolumnname in self.__data[i]:
                self.__data[i][self.__timecolumnname] /= self.__timerescale
        
        # future pickle implementations might use that
        self.__store_kwargs = kwargs


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
            if not wl[1] in idm:
                idm[wl[1]] = list()
            idm[wl[1]].append(i)
        for keys in idm.keys():
            idm[keys] = np.array(idm[keys],dtype=int)
        return idm
    

    def ExcludeLabel(self,label):
        self.__excludelabels.append(label)

    def SetMaxTime(self,maxtime):
        self.__maxtime = maxtime

    def InstantGrowthRates(self,dropID,signalcolumn,logsignal = True):
        growthrates = list()
        if self.__timecolumnname in self.__data[dropID] and signalcolumn in self.__data[dropID]:
            if self.__splitBackForthTrajectories:
                t0 = self.__data[dropID][self.__timecolumnname][0::2]
                t1 = self.__data[dropID][self.__timecolumnname][1::2]
                s0 = self.__data[dropID][signalcolumn][0::2]
                s1 = self.__data[dropID][signalcolumn][1::2]
                
                if not self.__maxtime is None:
                    s0 = s0[t0 < self.__maxtime]
                    s1 = s1[t1 < self.__maxtime]
                    t0 = t0[t0 < self.__maxtime]
                    t1 = t1[t1 < self.__maxtime]
                
                if logsignal:
                    t0 = t0[s0>0]
                    t1 = t1[s1>0]
                    s0 = np.log(s0[s0>0])
                    s1 = np.log(s1[s1>0])
                
                if len(t0) > self.__fitpoints:
                    gr0  = np.array([[LMSQ(t0[i:i+self.__fitpoints],s0[i:i+self.__fitpoints])[0][1],np.mean(t0[i:i+self.__fitpoints])] for i in range(len(t0)-self.__fitpoints)],dtype= np.float)
                else:
                    gr0 = None
                
                if len(t1) > self.__fitpoints:
                    gr1 = np.array([[LMSQ(t1[i:i+self.__fitpoints],s1[i:i+self.__fitpoints])[0][1],np.mean(t0[i:i+self.__fitpoints])] for i in range(len(t1)-self.__fitpoints)],dtype= np.float)
                else:
                    gr1 = None
                    
                ret = (gr0,gr1)
            else:
                t0 = self.__data[i][self.__timecolumnname]
                s0 = self.__data[i][signalcolumn]
                
                if not self.__maxtime is None:
                    s0 = s0[t0 < self.__maxtime]
                    t0 = t0[t0 < self.__maxtime]
                
                if logsignal:
                    t0 = t0[s0>0]
                    s0 = np.log(s0[s0>0])
                
                if len(t0) > self.__fitpoints:
                    gr0 = np.array([[LMSQ(t0[i:i+self.__fitpoints],s0[i:i+self.__fitpoints])[0][1],np.mean(t0[i:i+self.__fitpoints])] for i in range(len(t0)-self.__fitpoints)],dtype= np.float)
                else:
                    gr0 = None
                ret = gr0
        return ret
                

    def GrowthRatesList(self,label,signalcolumn,maxfrac = .95,logsignal = True):
        growthrates = list()
        if label in self.labels:
            #fp = open('growthrates-{}'.format(label),'w')
            for dropID in self.__idropmap[label]:
                if self.__splitBackForthTrajectories:
                    gr0,gr1 = self.InstantGrowthRates(dropID,signalcolumn,logsignal)
                    if not gr0 is None:
                        mgr0 = MeanUpperFrac(gr0[:,0])
                        if not mgr0 is None:
                            growthrates.append(mgr0)
                    if not gr1 is None:
                        mgr1 = MeanUpperFrac(gr1[:,0])
                        if not mgr1 is None:
                            growthrates.append(mgr1)
                else:
                    gr0 = self.InstantGrowthRates(dropID,signalcolumn,logsignal)
                    if not gr0 is None:
                        mgr0 = MeanUpperFrac(gr0[:,0])
                        if not mgr0 is None:
                            growthrates.append(mgr0)
        
        growthrates = np.array(growthrates,dtype=np.float)
        growthrates = growthrates[~np.isnan(growthrates)]
        return growthrates

    def DropID_From_Label(self,label):
        if label in self.labels:
            return self.__idropmap[label]
        else:
            return None

    def data(self,dropID):
        return self.__data[dropID]

    def __getattr__(self,key):
        if key == 'labels':
            return np.sort([label for label in self.__idropmap.keys() if label.upper() not in self.__excludelabels])
        elif key == 'columns':
            return self.__columns
        #else:
            #super(DropletDataNew,self).__getattr__(key)
            
    
    def __getitem__(self,key):
        if key in self.__idropmap.keys():
            for dropID in self.__idropmap[key]:
                yield self.__data[dropID]
        else:
            raise KeyError('"{}" not in allowed labels: ('.format(key) + ', '.join(self.__idropmap.keys()) + ')')



def LMSQ(x,y):
    # least mean squares estimator
    # for a linear interpolation, including covariance matrix for the estimated parameters
    # A = ( 1  ... 1  )
    #     ( x1 ... xn )
    #
    # y = ( y1 ... yn ).T
    #
    # p = (a b).T
    #
    # E[p]     = inv(A.T * A) * A.T * y
    # Cov[p,p] = sigma2 * inv(A.T * A)
    # sigma2   = E[ ( y - A*E[p] )^2 ]
    n   = len(x)
    sx  = np.sum(x)
    sy  = np.sum(y)
    sxx = np.dot(x,x)
    sxy = np.dot(x,y)
    syy = np.dot(y,y)
    
    # estimate parameters
    denom  = (n*sxx-sx*sx)
    b      = (n*sxy - sx*sy)/denom
    a      = (sy-b*sx)/n
    estim  = np.array([a,b],dtype=np.float)

    # estimate covariance matrix of estimated parameters
    sigma2 = syy + n*a*a + b*b*sxx + 2*a*b*sx - 2*a*sy - 2*b*sxy # variance of deviations from linear line
    cov    = sigma2 / denom * np.array([[sxx,-sx],[-sx,n]],dtype=np.float)

    return estim,cov


def MeanUpperFrac(x,upperfrac = .95,returnnumber = False):
    if not x is None:
        # delete all NaN values
        xx = np.array([a for a in x if not np.isnan(a)])
        if len(xx) > 0:
            uf_val = np.max(xx) * upperfrac
            m = np.mean(xx[xx >= uf_val])
            n = len(xx[xx >= uf_val])
            if returnnumber:
                return m,n
            else:
                return m
        else:
            return None
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


