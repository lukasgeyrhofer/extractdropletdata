#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys,math
import argparse
import os

import millidrop_dataclass_new as mdc

signaldata = {'WT':'fluo_2_mean', 'PVDS': 'fluo_3_median' }


parser = argparse.ArgumentParser()
parser.add_argument("-i","--rootdir")
parser.add_argument("-E","--excludelabels",nargs="*",default = ['Empty'])
args = parser.parse_args()

dirlist = np.sort([d for d in os.listdir(args.rootdir) if os.path.isdir(d)])
excludelabels = [l.upper() for l in args.excludelabels]

for expdir in dirlist:
    print expdir

    templatefile = os.path.join(expdir,'analysis/template.csv')
    dropdir      = os.path.join(expdir,'analysis/droplets/')

    kwargs = {  'templatefile':templatefile,
                'infiles':dropdir,
                'excludelabels':excludelabels}

    data = mdc.DropletDataNew(**kwargs)

    curstrain = None
    for a in signaldata.keys():
        if a in expdir.upper():
            curstrain = a
    if not curstrain is None:
        for label in data.labels:
            if signaldata[curstrain] in data.columns:
                gr = data.GrowthRatesList(label,signaldata[curstrain])
                print '   {:20s} {:.6f} {:.6f} {:.6f} {:4d}'.format(label,np.mean(gr),np.std(gr),np.median(gr),len(gr))



