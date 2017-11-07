#!/usr/bin/env python

import numpy as np
import argparse
import sys,math

parser = argparse.ArgumentParser()
parser.add_argument("-t","--templatefile",default=None)
parser.add_argument("-o","--outfile",default=None)
parser.add_argument("-d","--delimiter",default=',')
args = parser.parse_args()


try:
    fp = open(args.templatefile)
except:
    raise IOError("could not open templatefile '{}'.".format(args.templatefile))

if args.outfile is None:
    out = sys.stdout
else:
    try:
        out = open(args.outfile,"w")
    except:
        raise IOError("could not open file '{}' for writing.".format(args.outfile))

first = True
for line in fp.readlines():
    if len(line) > 1:
        if line[0] == "#":
            out.write(line)
        else:
            values = line.split(args.delimiter)
            if first:
                IDwell        = [x.lower() for x in values].index("well")
                IDdescription = [x.lower() for x in values].index("description")
                first = False
            else:
                values[IDdescription] += "-{}".format(values[IDwell])
            out.write(args.delimiter.join(values))
        
out.close()



