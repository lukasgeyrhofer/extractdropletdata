#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys,math
import argparse
import os

import millidrop_dataclass_new as mdc


def writetofile(fp,x):
    try:
        for a,b in x:
            fp.write('{:.6f} {:.6f}\n'.format(b,a))
        fp.write('\n')
    except:
        pass


signaldata = {'WT':'fluo_2_median', 'PVDS': 'fluo_3_median', 'ECOLI': 'fluo_2_mean'}


parser = argparse.ArgumentParser()
parser.add_argument("-i","--rootdir")
parser.add_argument("-E","--excludelabels",nargs="*",default = ['Empty','KB','CAA','LB'])
parser.add_argument("-m","--maxtime",type=float,default=None)
args = parser.parse_args()

dirlist = np.sort([d for d in os.listdir(args.rootdir) if os.path.isdir(d) and d[0] != '.'])
excludelabels = [l.upper() for l in args.excludelabels]

alldata = list()
labelList = list()

for expdir in dirlist:
    print expdir

    templatefile = os.path.join(expdir,'analysis/template.csv')
    dropdir      = os.path.join(expdir,'analysis/droplets/')

    kwargs = {  'templatefile':templatefile,
                'infiles':dropdir,
                'excludelabels':excludelabels}
    
    data = mdc.DropletDataNew(**kwargs)
    if not args.maxtime is None:
        data.SetMaxTime(args.maxtime)
        
    curstrain = None
    for a in signaldata.keys():
        if a in expdir.upper():
            curstrain = a
    if not curstrain is None:
        for label in data.labels:
            if signaldata[curstrain] in data.columns:
                fp = open(os.path.join(expdir,'growthrates-{}.txt'.format(label)),'w')
                for dropID in data.DropID_From_Label(label):
                    gr0,gr1 = data.InstantGrowthRates(dropID,signaldata[curstrain])
                    writetofile(fp,gr0)
                    writetofile(fp,gr1)
                fp.close()
                gr = data.GrowthRatesList(label,signaldata[curstrain])
                if not gr is None:
                    if len(gr) > 0:
                        print '   {:20s} {:10s} {:.6f} {:.6f} {:.6f} {:4d}'.format(label,signaldata[curstrain],np.mean(gr),np.std(gr),np.median(gr),len(gr))


data=[]
labelList=[]
for path in folder: #[folder[0]]:
    a=path.split('/')
    pathGrate=path+'growthDerivee/'
    tmp=pd.read_csv(pathGrate+'gRate_nbpt5')

    for experiment in tmp.columns:
        if experiment != 'Unnamed: 0':
            x=np.array(tmp[experiment])
            x=x[~np.isnan(x)]
            data.append(x)
            labelList.append(a[-3]+experiment)
