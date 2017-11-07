#!/usr/bin/env python3

import numpy as np
import argparse
import sys,math
import os
import pandas as pd
from itertools import groupby


class DropletData(object):
    def __init__(self, **kwargs):

        self.__restrictedvaluesforparameters = list([None,""])
        
        # list of all filenames, usually in droplets folder from MilliDrop pipeline, ie. droplets/0000.csv, droplets/0001.csv, ...
        self.__infiles                    = kwargs.get("infiles",None)
        # templatefile to assign the labels from different wells
        self.__templatefile               = kwargs.get("templatefile",None)
        # if restrictions are loaded from an external file
        self.__restrictionfile            = kwargs.get("restrictionfile",None)

        # which channels should be loaded, by default load everything
        self.__datacolumns                = kwargs.get("datacolumns",None)

        # downstream processing, not actively accessed by this data object
        # used in several of the derivative scripts that use this object, only used to store this value
        self.__timerescale                = kwargs.get("timerescale",3.6e3)
        # store folder to output files
        self.__outbasename                = kwargs.get("outbasename","")
        if self.__outbasename is None:
            self.__outbasename = ""
        
        # further options when loading the data
        self.__splitBackForthTrajectories = kwargs.get("SplitBackForthTrajectories",True)
        self.__snakelikeloading           = kwargs.get("SnakeLikeLoading",True)
        self.__nonhiccuploading           = kwargs.get("NonHiccupLoading",False)
        self.__ignoreadditionaldroplets   = kwargs.get("IgnoreAdditionalDroplets",False)
                                                       

        if self.__infiles is None:
            raise IOError("datafiles required")
        if len(self.__infiles) < 1:
            raise IOError("datafiles required")
        
        
        if not self.__templatefile is None:
            self.load_templatefile(self.__templatefile)
        else:
            # if no templatefile, group everything under label 'default' with 'unkown' well ID
            self.__droplettype = np.repeat("default",len(infiles))
            self.__well        = np.repeat("unknown",len(infiles))
        
        # ===============================================================
        # = variable initialization
        # ===============================================================
        self.__listoftypes         = list(set(self.__droplettype))  # list of all possible "experiments", ie. labels of different wells
        self.__data                = dict()                         # store all trajectories, each keys of the dictionary are the experiment labels, "values" of such an entry is list of np-arrays
        self.__welldata            = dict()
        self.__datarestrictions    = list()                         # list of all restrictions
        self.__permittedoperations = {"min":2,"max":2,"end":2,"start":2,"exclude":1}
                                                                    # hardcoded allowed restriction types: min/max chop trajectories if the fall below or rise above the threshold value, end/start is for time


        # ===============================================================
        # = iterate over all separate droplet files
        # ===============================================================
        for filename in self.__infiles:
            dropletID = self.filename_to_dropletID(filename)
            try:
                tmpdata = np.genfromtxt(filename,names = True, delimiter = ',', dtype = float)
            except:
                print("Error while loading file '{:s}'. Continuing ...".format(filename))
                continue
            
            if self.__datacolumns is None:
                self.__datacolumns = tmpdata.dtype.names
            
            self.add_trajectory(dropletID, tmpdata)
        
        
        # ===============================================================
        # = restrictions on data
        # ===============================================================
        
        if not self.__restrictionfile is None:
            self.load_restrictions_from_file(self.__restrictionfile)
        


    # ===============================================================
    # = generate list of all types of experiments from templatefile =
    # ===============================================================
    def load_templatefile(self,templatefile = None, templatefileseparator = ','):
        try:
            templatedata = pd.read_csv(templatefile, comment = "#", sep = templatefileseparator)
        except:
            raise IOError("could not open templatefile '{}'".format(templatefile))
        
        if 'well' in templatedata and 'description' in templatedata and 'droplet_number' in templatedata:
            well           = np.array(templatedata['well'],dtype=str)
            description    = np.array(templatedata['description'],dtype=str)
            droplet_number = np.array(templatedata['droplet_number'],dtype=int)
            
        else:
            raise ValueError("relevant data not found in templatefile. need columns 'well', 'description' and 'droplet_number'")
        
        
        # generate correct order of loading from microwell plate
        if 'order' in templatedata:
            # newer versions of templatefile seem to have an 'order' column, ignore flag 'snakelikeloading' if this is present
            order = np.array(templatedata['order'],dtype=int)
            order, well, description, droplet_number = zip(*sorted(zip(order,well,description,droplet_number))) # sort all columns with respect to 'order'

        elif self.__snakelikeloading:
            # construct same order list as above by hand
            order = None
            totalcount = 0
            all_rows = [x[0] for x in well]
            all_row_separate = [list(j) for i,j in groupby(all_rows)]
            for a in all_row_separate:
                direction   = 2 * (ord(a[0])%2) - 1
                countrow    = len(a)
                order       = self.concat(order, np.arange(totalcount,totalcount + countrow)[::direction])
                totalcount += countrow
            order, well, description, droplet_number = zip(*sorted(zip(order,well,description,droplet_number))) # sort all columns with respect to 'order'
        else:
            # otherwise just use default order
            pass
        
        
        # use this correct order to generate lists with experiments and wells, which will be later matched to a dropletID
        self.__droplettype = None
        self.__well        = None
        if not self.__nonhiccuploading:
            for i in np.arange(start = 0,stop = len(well),step = 2):
                assert droplet_number[i] == droplet_number[i+1]
                self.__droplettype = self.concat(self.__droplettype, np.repeat([[description[i],description[i+1]]],droplet_number[i],axis=0).flatten())
                self.__well        = self.concat(self.__well,        np.repeat([[well[i],well[i+1]]],droplet_number[i],axis=0).flatten())
        else:
            for w,d,n in zip(well,description,droplet_number):
                self.__droplettype = self.concat(self.__droplettype, np.repeat(d,n))
                self.__well        = self.concat(self.__well,        np.repeat(w,n))
                                                                   

        
    
        

    # ===============================================================
    # = load data from a single dropletfile and add to internal datastructure
    # ===============================================================
    def add_trajectory(self,dropletID = None, data = None):
        if data is None:
            raise ValueError("need timestamps for experimental data")
        
        

        for column in self.__datacolumns:
            if not column in data.dtype.names:
                raise ValueError("DropletID: {} --- No field of name '{}'. Possible values are ('".format(dropletID,column) + "', '".join(data.dtype.names) + "')")
        try:
            dropletLabel = self.dropletID_to_label(dropletID)
            dropletWell  = self.dropletID_to_well(dropletID)
        except:
            if not self.__ignoreadditionaldroplets:
                raise KeyError("Error while assigning Label and/or Well to droplet")
            else:
                return
        
        if not dropletLabel in self.__data:
            self.__data[dropletLabel]     = list()
            self.__welldata[dropletLabel] = list()
        
        if self.__splitBackForthTrajectories:
            newdatablock0 = dict()
            newdatablock1 = dict()
            for column in self.__datacolumns:
                newdatablock0[column] = np.array(data[column][0::2])
                newdatablock1[column] = np.array(data[column][1::2])                                                    
            self.__data[dropletLabel].append(newdatablock0)
            self.__data[dropletLabel].append(newdatablock1)
            self.__welldata[dropletLabel].append(dropletWell)
            self.__welldata[dropletLabel].append(dropletWell)
        else:
            newdatablock0 = dict()
            for column in self.__datacolumns:
                newdatablock0[column] = np.array(data[column])
            self.__data[dropletLabel].append(newdatablock0)
            self.__welldata[dropletLabel].append(dropletWell)
    
    
    
    # ===============================================================
    # = helper routines
    # ===============================================================
    
    def concat(self,list1 = None,list2 = None, direction1 = 1, direction2 = 1):
        if (list1 is None) and (list2 is None):
            return None
        elif (list1 is None):
            return list2[::direction2]
        elif (list2 is None):
            return list1[::direction1]
        else:
            return np.concatenate([list1[::direction1],list2[::direction2]])


    def filename_to_dropletID(self,filename):
        # in current implementation, filename is just 'dropletnumber.CSV'
        # change here if necessary
        r = int(os.path.basename(filename).split('.')[0])
        return r


    def dropletID_to_label(self,dropletID = None):
        if dropletID is None:
            raise KeyError("dropletID is None")
        try:
            return self.__droplettype[dropletID]
        except:
            raise KeyError("did not find experiment type for droplet '{}'".format(dropletID))
    
    
    def dropletID_to_well(self,dropletID = None):
        if dropletID is None:
            raise KeyError("dropletID is None")
        try:
            return self.__well[dropletID]
        except:
            raise KeyError("did not find well type for droplet '{}'".format(dropletID))
    
    # ===============================================================
    # = output of specific information from internal datastructure
    # ===============================================================

    def CountTrajectories(self):
        r = dict()
        for key in self.__listoftypes:
            r[key] = len(self.__data[key])
        return r
    
    
    def getDropletTypes(self):
        return self.__droplettype
    
    
    def getChannelIndex(self,key):
        if key in self.__datacolumns:
            return self.__datacolumns.index(key)
        else:
            return None
    
    
    # ===============================================================
    # routines to work with restricted data (ie. a lower cutoff for the droplet signal, or a maximal time)
    # ===============================================================
    
    def set_restriction(self,restrictiontype, applies_to = None, params = None):
        if restrictiontype in self.__permittedoperations.keys():
            if applies_to is None:
                applies_to = "all"
            if (applies_to in self.__droplettype) or (applies_to in set(self.__well)) or (applies_to == "all"):
                if len(params) == self.__permittedoperations[restrictiontype]:
                    newrestriction = list([restrictiontype,applies_to,params])
                    self.__datarestrictions.append(newrestriction)
        else:
            raise ValueError("cannot apply restriction to data (column: '{:s}', operation: {:s})".format(column,operation))


    def remove_all_restrictions(self):
        self.__datarestrictions = list()
        
        
    def load_restrictions_from_file(self,filename):
        try:
            fp = open(filename,"r")
        except:
            raise IOError("could not open file '{:s}' to load restrictions".format(filename))
        for line in fp.readlines():
            values = line.split()
            if len(values) >= 2:
                if values[0] in self.__permittedoperations.keys():
                    if (values[1] in self.__droplettype) or (values[1] in set(self.__well)) or (values[1] == "all"):
                        if len(values) == 2 + self.__permittedoperations[values[0]]:
                            self.set_restriction(values[0],values[1],values[2:])
        fp.close()
    
    
    def write_restrictions_to_file(self,filename = None):
        try:
            if filename is None:
                fp = sys.stdout
            else:
                fp = open(filename,"w")
        except:
            raise IOError("could not open file '{:s}' to save restrictions".format(filename))
        for restriction in self.__datarestrictions:
            fp.write(restriction[0] + " " + restriction[1] + " " + " ".join([str(x) for x in restriction[2]]))
        if not filename is None:
            fp.close()

    # ===============================================================
    # all output of data is funneled through this routine
    # only return datapoints that match all criteria in the restrictionfile
    # ===============================================================
    
    def pattern_and(self,pattern1,pattern2):
        if pattern1 is None:
            return pattern2
        else:
            return np.logical_and(pattern1,pattern2)
    
    def restricted_data(self,key):
        r = list()
        for i in range(len(self.__data[key])):
            datablock = self.__data[key][i]
            well      = self.__welldata[key][i]
            rdata = None
            for column in self.__datacolumns:
                if rdata is None:
                    rdata = np.array([datablock[column]])
                else:
                    rdata = np.concatenate([rdata,np.array([datablock[column]])],axis = 0)
            rdata = np.transpose(rdata)
            pattern = None
            if len(self.__datarestrictions) > 0:
                for restriction in self.__datarestrictions:
                    if (restriction[1] == key) or (restriction[1] == well) or (restriction[1].lower() == "all"):
                        if   str(restriction[0]).lower() == "max":
                            try:    pattern = self.pattern_and(pattern, datablock[restriction[2][0]] < float(restriction[2][1]) )
                            except: continue
                        elif str(restriction[0]).lower() == "min":
                            try:    pattern = self.pattern_and(pattern, datablock[restriction[2][0]] > float(restriction[2][1]) )
                            except: continue
                        elif str(restriction[0]).lower() == "end":
                            try:    pattern = self.pattern_and(pattern, datablock[restriction[2][0]] < (1-restriction[2][1]) * datablock[restriction[2][0]][-1] )
                            except: continue
                        elif str(restriction[0]).lower() == "start":
                            try:    pattern = self.pattern_and(pattern, datablock[restriction[2][0]] > (1+restriction[2][1]) * datablock[restriction[2][0]][0] )
                            except: continue
                        elif str(restriction[0]).lower() == "exclude":
                            if ((str(restriction[2]).lower() == "experiment") and (key == str(restriction[1]))) or ((str(restriction[2]).lower() == "well") and (well == str(restriction[1]))):
                                shape = np.shape(rdata)
                                pattern = np.repeat(np.repeat([[False]],shape[0],axis=0),shape[1],axis=1)
                            else:
                                continue

                pattern = np.transpose(np.repeat([pattern],len(rdata[0]),axis = 0))
                datapattern = rdata[pattern]
                if len(datapattern) > 0:
                    rdata   = np.reshape(datapattern,(len(datapattern)/len(self.__datacolumns),len(self.__datacolumns)))
                else:
                    rdata   = datapattern
            if len(rdata) > 0:
                r.append(rdata)
        return r
          
    
    # ===============================================================
    # = object logic                                                =
    # ===============================================================
    
    # after object initialization, this is how the dataobject reacts
    # example code lines start with '>'
    # > data = DropletData()
    
    # get all trajectories of a certain experiment as list,
    # when EXPERIMENTLABEL is a string also contained in TEMPLATEFILE as label for a well
    # > data[EXPERIMENTLABEL]
    def __getitem__(self,key):
        if key in self.__listoftypes:
            return self.restricted_data(key)
    
    # iterate over the whole dataset
    # > for ExperimentLabel, Trajectories in data:
    # >     print ExperimentLabel
    # >     for trajectory in Trajectories:
    # >         print trajectory
    def __iter__(self):
        for key in self.__listoftypes:
            data = self[key]
            if len(data) > 0:
                yield key,data

    # catch any specific attributes of the DropletData object
    # so far only single attribute implemented
    # > print data.datacolumns
    def __getattr__(self,key):
        if   key == "datacolumns":
            return self.__datacolumns
        elif key == "outbasename":
            return self.__outbasename
        elif key == "listoftypes":
            return self.__listoftypes
        else:
            super(DropletData,self).__getattr__(key)

# helper script to add all cmdline parameters
def AddCommandLineParameters(parser):
    ioparser = parser.add_argument_group(description = "==== I/O parameters ====")
    
    ioparser.add_argument("-i", "--infiles",         nargs="*")
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




# working example of loading files and printing them again
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--infiles",nargs="*")
    parser.add_argument("-T","--templatefile",default=None)
    parser.add_argument("-S","--splitBackForthTrajectories",default=False,action="store_true")
    args = parser.parse_args()
    
    # define data object
    data = DropletData(infiles = args.infiles, templatefile = args.templatefile, splitBackForthTrajectories = args.splitBackForthTrajectories)

    # in general, data is obtained by iterating over sets of experimental labels with their trajectories
    for dropletLabel, Trajectories in data:
        print("'{:s}' with {:d} trajectories".format(dropletLabel, len(Trajectories)))
        for trajectory in Trajectories:
                print(trajectory)
                
    # should also work like that:
    # given that the label 'SBW25mC-10' is defined in the templatefile
    # this is the same as the 'Trajectories' in the loop above
    #if 'SBW25mC-10' in data.listoftypes:
        #print data['SBW25mC-10']


# run through working example, if program is executed directly
# don't run, if just loaded for class definitions
if __name__ == "__main__":
    main()






