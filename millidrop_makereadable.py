#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import argparse
import sys,math


import millidrop_dataclass as mdc

parser = argparse.ArgumentParser()
parser = mdc.AddCommandLineParameters(parser)

dparser = parser.add_argument_group(description = "==== Data parameters ====")
dparser.add_argument("-c","--signalcolumn",default = 0,type = float)
dparser.add_argument("-x","--lowertimecutoff",default = None,type = float)
dparser.add_argument("-X","--uppertimecutoff",default = None,type = float)
dparser.add_argument("-s","--lowersignalcutoff",default = None,type = float)
dparser.add_argument("-S","--uppersignalcutoff",default = None,type = float)
args = parser.parse_args()

data = mdc.DropletData(**vars(args))

if not args.lowertimecutoff is None:
    data.set_restriction("time","min",args.lowertimecutoff)
if not args.uppertimecutoff is None:
    data.set_restriction("time","max",args.uppertimecutoff)
if not args.lowersignalcutoff is None:
    data.set_restriction(data.datacolumns[args.signalcolumn],"min",args.lowersignalcutoff)
if not args.uppersignalcutoff is None:
    data.set_restriction(data.datacolumns[args.signalcolumn],"max",args.uppersignalcutoff)


for label,trajectories in data:
    i = 0
    for traj in trajectories:
        if len(traj) > 0:
            filename = data.outbasename + "{:s}-{:04d}.data".format(label,i)
            np.savetxt(filename,traj,fmt = '%.6e', delimiter = ' ')
            i += 1
    print("{:12s}: saving {:4d} trajectories to readable format".format(label,i))
