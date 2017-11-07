#!/usr/bin/env python

import numpy as np
import argparse
import sys,math


def FirstIndex(label,data):
    r = [x[0] for x in data]
    try:    index = r.index(label)
    except: index = None
    return  index


def ReverseData(data):
    r = list()
    for x in data:
        if ord(x[0])%2 == 1:
            r.append(x)
        else:
            index = FirstIndex(x[0],r)
            if index is None:   r.append(x)
            else:               r.insert(index,x)
    return r


def LoadData(filename):
    try:    fp = open(filename,"r")
    except: raise IOError("could not open template '{:s}'".format(filename))
    names = None
    data  = list()
    for line in fp.readlines():
        if (line.strip())[0] != "#":
            values = line.split(',')
            if names is None:
                names = [v.strip() for v in values]
                try:
                    wellID = names.index('well')
                    descriptionID = names.index('description')
                    droplet_numberID = names.index('droplet_number')
                except:
                    raise ValueError("template file does not contain necessary data")
            else:
                data.append([values[wellID][0],values[wellID][1:],values[descriptionID],int(values[droplet_numberID])])
    return names,data


def WriteData(filename,data,names,delimiter = ','):
    try:    fp = open(filename,"w")
    except: raise IOError("could not open file '{}' for writing".format(filename))
    print >> fp,delimiter.join(names)
    for x in data:
        print >> fp,"{}{},{},{}".format(*x)
    fp.close()


def ReplaceEmptyWells(data,emptylabel):
    r = list()
    for x in data:
        if x[2] != emptylabel:
            r.append(x)
        else:
            for i in range(x[3]):
                r.append([x[0],x[1],emptylabel,1])
    return r


def ComputeDistancesToNonEmpty(data,index,emptylabel):
    i = index
    while data[i][2] == emptylabel: i -= 1
    if i >= 0:          prevNonEmpty = (index - i,data[i][2])
    else:               prevNonEmpty = (None,None)
        
    i = index
    while data[i][2] == emptylabel: i += 1
    if i < len(data):   nextNonEmpty = (i - index,data[i][2])
    else:               nextNonEmpty = (None,None)
            
    return prevNonEmpty,nextNonEmpty


def WriteDistanceToLabel(data,emptylabel):
    r = list()
    for i in range(len(data)):
        if data[i][2] != emptylabel:
            r.append(data[i])
        else:
            ne,la = ComputeDistancesToNonEmpty(data,i,emptylabel)
            if ne[0] is None:
                r.append([data[i][0],data[i][1],data[i][2] + '-{:s}-{:03d}'.format(la[1],la[0]),data[i][3]])
            elif la[0] is None:
                r.append([data[i][0],data[i][1],data[i][2] + '-{:s}-{:03d}'.format(ne[1],ne[0]),data[i][3]])
            elif la[0] <= ne[0]:
                r.append([data[i][0],data[i][1],data[i][2] + '-{:s}-{:03d}'.format(la[1],la[0]),data[i][3]])
            elif la[0] > ne[0]:
                r.append([data[i][0],data[i][1],data[i][2] + '-{:s}-{:03d}'.format(ne[1],ne[0]),data[i][3]])
    return r


def main():    
    parser = argparse.ArgumentParser()
    parser.add_argument("-t","--templatefile",default=None,required=True)
    parser.add_argument("-e","--emptylabel",default="KBIPTG",type=str)
    parser.add_argument("-S","--snakelikeloading",default=True,action="store_false")
    parser.add_argument("-o","--outfile",default=None,required=True)
    args = parser.parse_args()

    names,data = LoadData(args.templatefile)

    if args.snakelikeloading:
        data = ReverseData(data)
    data = ReplaceEmptyWells(data,args.emptylabel)
    data = WriteDistanceToLabel(data,args.emptylabel)
    if args.snakelikeloading:
        data = ReverseData(data)

    WriteData(args.outfile,data,names)



if __name__ == "__main__":
    main()
