#!/usr/bin/env python

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
        
        # further options when loading the data
        self.__splitBackForthTrajectories = kwargs.get("SplitBackForthTrajectories",True)
        
        
        self.__dropmap  = self.LoadDropMap(self.__templatefile)
        self.__idropmap = self.InvertDropmap(self.__dropmap)
        
        self.__data     = dict()
        for i,fn in enumerate([self.__dropdir + f for f in os.listdir(self.__dropdir) if os.path.isfile(os.path.join(self.__dropdir,f))]):
            self.__data[i] = pd.read_csv(fn)

    def LoadDropMap(self,filename):
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
            if not idm.has_key(wl[1]):
                idm[wl[1]] = list()
            idm[wl[1]].append(i)
        for keys in idm.iterkeys():
            idm[keys] = np.array(idm[keys],dtype=int)
        return idm
