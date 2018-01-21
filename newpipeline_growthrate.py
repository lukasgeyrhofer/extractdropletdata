#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys,math
import argparse
import os

import millidrop_dataclass_new as mdc

signaldata = {'WT':'fluo_2_median', 'PVDS': 'fluo_3_median' }


parser = argparse.ArgumentParser()
parser.add_argument("-i","--rootdir")
parser.add_argument("-E","--discardlabels",nargs="*",default = [])
args = parser.parse_args()

dirlist = np.sort([d for d in os.listdir(args.rootdir) if os.path.isdir(d)])
discardlabels = [l.upper() for l in args.discardlabels]


for expdir in dirlist:
    print expdir

    templatefile = os.path.join(expdir,'analysis/template.csv')
    dropdir      = os.path.join(expdir,'analysis/droplets/')

    kwargs = {  'templatefile':templatefile,
                'infiles':dropdir}

    data = mdc.DropletDataNew(**kwargs)

    for a in signaldata.keys()

    for label in [l for l in data.labels if not l.upper() in discardlabels]:
        for trajectory in data[label]:
            columns = [str(a) for a in trajectory.keys()]

            

            time = trajectory['time'] / 3.6e3
            signal = np.log(trajectory['test'])
    














