#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==================================================================== #
#    Extract data from droplet train                                   #
# ==================================================================== #
#                                                                      #
#    Takes time series data and extracts info for each droplet         #
#    by finding increases in signal that corresponds to a              #
#    droplet, then fitting fit fuctions through the time series        #
#    to get properties of the droplet.                                 #
#                                                                      #
#    Currently implemented fit functions:                              #
#      * rectangular                                                   #
#      * Gaussian                                                      #
#      * sigmoid                                                       #
#    Each fit yields                                                   #
#      * position of droplet center in time                            #
#      * maximum value of signal                                       #
#      * base value of signal                                          #
#      * parameter related to width of droplet                         #
#      * (for sigmoid only:) steepness of signal at flanks             #
#                                                                      #
#    Usage:                                                            #
#      ./extract_data.py -i DATAFILE -m 0.15 -M 0.3                    #
#                                                                      #
#    Use option --help to see improved description                     #
#                                                                      #
#                                                                      #
#    Lukas Geyrhofer, 2016                                             #
#                                                                      #
#                                                                      #
# ==================================================================== #



import numpy as np
import argparse
import sys,math
from scipy.optimize import curve_fit


# ==================================================================== #
# data structure for droplets                                          #
# ==================================================================== #
#   contains fitdata, iterating over class yields the data from        #
#   all fits currently stored in the object                            #
# ==================================================================== #
class droplet:
    def __init__(self):
        self.__dropletcount = 0
        
        self.__dropletmax  = None
        self.__dropletbase = None
        self.__droplettime = None
        self.__dropletdataindex = None
        
        self.__fitdata = {}
        self.__fitkeys = []

    
    def add_droplet(self,maxval,baseval,time,idx = 0):
        if self.__dropletcount > 0:
            self.__dropletmax  = np.concatenate((self.__dropletmax ,np.array([maxval])  ))
            self.__dropletbase = np.concatenate((self.__dropletbase,np.array([baseval]) ))
            self.__droplettime = np.concatenate((self.__droplettime,np.array([time])    ))
            self.__dropletdataindex = np.concatenate((self.__dropletdataindex,np.array([idx]) ))
        else:
            self.__dropletmax  = np.array([maxval])
            self.__dropletbase = np.array([baseval])
            self.__droplettime = np.array([time])
            self.__dropletdataindex = np.array([idx])
        self.__dropletcount += 1
    
    
    def add_droplet_fitinfo(self,key,dropletID,values):
        c = len(values)
        if not key in self.__fitkeys:
            self.__fitdata[key] = np.zeros((self.__dropletcount,c))
            self.__fitkeys.append(key)
        if dropletID < self.__dropletcount:
            self.__fitdata[key][dropletID] = values
    
    # fitdata can be retrieved like a "dictionary"
    def __getitem__(self,key):
        if key in self.fitkeys:
            return self.__fitdata[key]
    
    def __getattr__(self,key):
        if key == "time":
            return self.__droplettime
        elif key == "max":
            return self.__dropletmax
        elif key == "base":
            return self.__dropletbase
        elif key == "index":
            return self.__dropletdataindex
        elif key == "count":
            return self.__dropletcount
        elif key == "data":
            return self.__droplettime,self.__dropletmax,self.__dropletbase
        elif key == "fitkeys":
            return self.__fitkeys
    
    # casting object as int returns number of droplets
    def __int__(self):
        return self.__dropletcount
            
    # iterating yields arrays with fitdata
    def __iter__(self):
        for key in self.__fitkeys:
            yield self[key]



# ==================================================================== #
# fitfunctions                                                         #
# ==================================================================== #
def gaussoffset(x,mean,stddev,sqrtmaxval,sqrtbaseval):
    return np.exp(-0.5*(x-mean)**2/stddev**2)*sqrtmaxval**2+sqrtbaseval**2
def stepoffset(x,t1,t2,sqrtmaxval,sqrtminval):
    r = np.ones(len(x))*sqrtminval**2
    i1 = ((x-t1)**2).argmin()
    i2 = ((x-t2)**2).argmin()
    r[i1:i2] = sqrtmaxval**2
    return r
def sigmoidoffset(x,mean,sqrtcrossover,sqrtsteepness,sqrtheight,sqrtbase):
    r = np.zeros(len(x))
    i1 = ((x-mean)**2).argmin()
    r[:i1] = sqrtheight**2 /(1.+np.exp(-(x[:i1]-mean+sqrtcrossover**2)/sqrtsteepness**2)) + sqrtbase**2
    r[i1:] = sqrtheight**2 /(1.+np.exp( (x[i1:]-mean-sqrtcrossover**2)/sqrtsteepness**2)) + sqrtbase**2
    return r


# ==================================================================== #
# class that reads time series and extracts data                       #
# ==================================================================== #
#   holds an object droplets (see class above) to store data           #
# ==================================================================== #
class timeseries:
    def __init__(self,filename,minmax_lowerthreshold = .15, minmax_upperthreshold = .3,maxfev=5000,verbose=False):
        self.__filename = filename
        try:
            self.__timeseriesdata = np.genfromtxt(filename)
        except:
            raise IOError
        
        self.__droplets = droplet()
        self.__datalength = len(self.__timeseriesdata[:,0])
        self.__minmax_lowerthreshold = minmax_lowerthreshold
        self.__minmax_upperhreshold = minmax_upperthreshold
        
        self.__maxfev = maxfev
    
    # extract maxima:
    #                                               *****
    #                          *****               *     *
    #                         *     *             *       *
    #                        *       *           *        *
    #    maxthreshold  ------O-------*-----------O---------*------------
    #                       *        *           *         *
    #                       *        *           *         *
    #                      *          *         *           *
    #                      *           *        *            *
    #    minthreshold  ----*-----------O-------*-------------O----------
    #                      *           *       *             *     *
    #                   ***             *******               *****
    #                        |         |         |           |
    #                        |inupper  |inlower  |inupper    |inlower
    #                        |         |         |           |
    #                ====== findmax ==>|<== findmax ========>|<== findmax =======
    #
    def get_droplet_data_minmax(self):
        inlower = True
        inupper = False
        enteredlower = 0
        for i in range(self.__datalength):
            if inlower and self.__timeseriesdata[i,1] > self.__minmax_upperhreshold:
                inupper = True
                inlower = False
                enteredupper = i
            if inupper and self.__timeseriesdata[i,1] < self.__minmax_lowerthreshold:
                inupper = False
                inlower = True
                
                baseidx = self.__timeseriesdata[enteredlower:i,1].argmin()+enteredlower
                baseval = self.__timeseriesdata[baseidx,1]
                maxidx  = self.__timeseriesdata[enteredlower:i,1].argmax()+enteredlower
                maxval  = self.__timeseriesdata[maxidx,1]
                time    = self.__timeseriesdata[maxidx,0]
                
                self.__droplets.add_droplet(maxval,baseval,time,maxidx)
                enteredlower = i
    
    def fit_step(self,startindex,finalindex):
        t = self.__timeseriesdata[startindex:finalindex,0]
        p = self.__timeseriesdata[startindex:finalindex,1]
        
        maxidx = p.argmax()
        minidx = p.argmin()
        
        p0 = np.array([t[maxidx]-0.03,t[maxidx]+0.03,np.sqrt(p[maxidx]),np.sqrt(p[minidx] if p[minidx]>0 else 0)])
        
        parameters,cov = curve_fit(stepoffset,t,p,p0=p0,maxfev=self.__maxfev)
        fitparam = np.array([0.5*(parameters[0]+parameters[1]),parameters[2]**2,parameters[3]**2,abs(parameters[1]-parameters[0])])
        return fitparam
        
    
    def fit_gaussian(self,startindex,finalindex):
        t = self.__timeseriesdata[startindex:finalindex,0]
        p = self.__timeseriesdata[startindex:finalindex,1]
        
        maxidx = p.argmax()
        minidx = p.argmin()
        
        p0 = np.array([t[maxidx],0.03,np.sqrt(p[maxidx]-p[minidx]),np.sqrt(p[minidx] if p[minidx]>0 else 0)])
        
        parameters,cov = curve_fit(gaussoffset,t,p,p0=p0,maxfev=self.__maxfev)
        fitparam = np.array([parameters[0],parameters[2]**2+parameters[3]**2,parameters[3]**2,abs(parameters[1])])
        return fitparam
    
    
    def fit_sigmoid(self,startindex,finalindex):
        t = self.__timeseriesdata[startindex:finalindex,0]
        p = self.__timeseriesdata[startindex:finalindex,1]
        
        maxidx = p.argmax()
        minidx = p.argmin()

        p0 = np.array([t[maxidx],np.sqrt(0.03),np.sqrt(0.003),np.sqrt(p[maxidx]-p[minidx]),np.sqrt(p[minidx] if p[minidx]>0 else 0)])
        
        parameters,cov = curve_fit(sigmoidoffset,t,p,p0=p0,maxfev=self.__maxfev)
        fitparam = np.array([parameters[0],parameters[3]**2+parameters[4]**2,parameters[4]**2,parameters[1]**2,parameters[2]**2])
        return fitparam
        
    
    def get_droplet_data_fit(self):
        for i in range(int(self.__droplets)):
            if i > 0:                       
                start = int((self.__droplets.index[i-1] + self.__droplets.index[i])/2)
            else:
                start = 0
            if i < int(self.__droplets) - 1:
                final = int((self.__droplets.index[i+1] + self.__droplets.index[i])/2)
            else:
                final = len(self.__timeseriesdata)
            
            valuesgauss = self.fit_gaussian(start,final)            
            self.__droplets.add_droplet_fitinfo("gauss",i,valuesgauss)
            
            #valuesstep = self.fit_step(start,final)
            #self.__droplets.add_droplet_fitinfo("step",i,valuesstep)
            
            #valuessigmoid = self.fit_sigmoid(start,final)
            #self.__droplets.add_droplet_fitinfo("sigmoid",i,valuessigmoid)

    def __str__(self):
        if int(self.__droplets) == 0:
            return "no data"
        else:
            dtime,dmax,dbase = self.__droplets.data
            
            s = ["%5d %e %e %e"%(i,dtime[i],dmax[i],dbase[i]) for i in range(int(self.__droplets))]
            for fitdata in iter(self.__droplets):
                for i in range(int(self.__droplets)):
                    for j in range(len(fitdata[i])):
                        s[i] += " %e"%fitdata[i,j]
            
            return "\n".join(s)
    
    def print_gnuplot_fits(self,filename = None):
        for fitkey in self.__droplets.fitkeys:
            if   fitkey == "gauss":     print "fitgauss(x,time,maxval,minval,width) = (maxval-minval)*exp(-(x-time)**2/(2*width**2))+minval"
            elif fitkey == "step":      print "fitstep(x,time,maxval,minval,width) = (x<=time-width/2?minval:(x<=time+width/2?maxval:minval))"
            elif fitkey == "sigmoid":   print "fitsigmoid(x,time,maxval,minval,width,steepness) = (x>=time?minval+(maxval-minval)/(1+exp((x-time-width)/steepness)):minval+(maxval-minval)/(1+exp((-x+time-width)/steepness)))"
            else:                       print "ERROR: fit not implemented"
        if filename != None:
            print "plot \"%s\" u 1:2 w lp lw 5 lc rgb \"#babdb6\", "%filename,
        else:
            print "plot ",
        print ", ".join([", ".join(["fit%s(x,"%fitkey + ",".join(["%e"%v for v in value]) + ")" for value in self.__droplets[fitkey]]) for fitkey in self.__droplets.fitkeys()])
        
        
            
    


# ==================================================================== #
#  main part of program                                                #
# ==================================================================== #
#   parse command line, instantiate timeseries class with              #
#   data file, get maxima, fit shape of data                           #
# ==================================================================== #
def main():
    # parse commandline options
    parser = argparse.ArgumentParser(description = "Extract data from time series of PMT measurements in droplets.")
    parser.add_argument("-i","--infile",help = "Data file containing two columns with time stamp and PMT measurement")
    parser.add_argument("-m","--minthreshold",type=float,default=.15 , help = "Used for finding maxima. Should be higher than all base values.")
    parser.add_argument("-M","--maxthreshold",type=float,default=.3, help = "Used for finding maxima. Should be smaller than all maximal values in droplets.")
    parser.add_argument("-v","--verbose",default=False,action="store_true",help = "Print plot statements for fit functions for each droplet in gnuplot")
    args = parser.parse_args()


    # read data
    data = timeseries(args.infile,args.minthreshold,args.maxthreshold)
    # iterate trough time series to get maxima (and thus roughly the time of each droplet)
    data.get_droplet_data_minmax()
    # use times from maxima to fit functions through the signal
    data.get_droplet_data_fit()

    # output
    if args.verbose:
        # full output to lot all fit functions for each droplet
        data.print_gnuplot_fits(filename = args.infile)
    else:
        # otherwise just print table of data
        print data


# ==================================================================== #
# python: program could be loaded as module, if executed as            #
# standalone, run main routine                                         #
# ==================================================================== #
if __name__ == "__main__":
    main()

