#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import sys,math

import millidrop_dataclass as mdc


parser = argparse.ArgumentParser()
parser = mdc.AddDropletParameters(parser)


args = parse.parse_args()



data = mdc.DropletData(**vars(args))


for label,trajectories in data:
    print label
