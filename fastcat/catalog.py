from __future__ import print_function, division

import numpy as np
import astropy.units as u
try:
    import h5py
except:
    print ("Cannot import h5py, if you try to write to h5 it won't work")

class Catalog(object):
    """ 
    Basic object to hold a catalog of observed astronomical objects.
    Intentially very simple for the time being.

    It holds a structured array which you can access directly.
    Eg. cat["ra"] will give you 1D array of ra coordinas. Valid names are:
    
    * "ra", "dec": float, ra,dec coordinates in radians
    * "z": float, redshift
    * "r": float, distance in Mpc/h
    * "rmag": float, rmagnitude
    * "e1","e2" : float, intrinsic ellipticity
    * "g1","g2" : float, shears
    It also has a placeholder for meta-data, which is empty at the moment

    On construction:
    Options
    -------
    N: int
       number of objects in the catalog
    meta: string
       string containing meta info
   """
    
    def __init__ (self, N, fields=['ra','dec','z'],dNdz=None, bz=None,window=None,
                  photoz=None,meta=None):
        self.data=np.zeros(N,dtype=map(lambda x:(x,np.float32),fields))
        self.dNdz=dNdz
        self.bz=bz
        self.window=window
        self.photoz=photoz
        self.meta=meta

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key,item):
        self.data[key]=item

    def readH5(self, fname):
        """ 
        Reads Catalog from H5 file, specified as argument
        """
        of=h5py.File(fname, "r")
        self.data=of["data"].value
        print("Fix reader")
        stop()
        print("Read %i entries from %s"%(len(self.data),fname))

    def dumpH5(self, fname):
        """ 
        Write Catalog to H5 file, specified as argument
        """
        of=h5py.File(fname, "w")
        if (self.meta):
            meta=of.create_dataset("meta",data=[])
            for v in self.meta.keys():
                meta.attrs[v]=self.meta[v]
        dset=of.create_dataset("objects", data=self.data, chunks=True,
                               shuffle=True,compression="gzip", compression_opts=9)
        if type(self.dNdz)!=type(None):
            dset=of.create_dataset("dNdz", data=self.dNdz)
        if type(self.bz)!=type(None):
            dset=of.create_dataset("bz", data=self.bz)
        if type(self.window)!=type(None):
            wfunc=of.create_dataset("window",data=[])
            for k,v in self.window.items():
                wfunc.attrs[k]=v
        if type(self.photoz)!=type(None):
            photoz=of.create_dataset("photoz",data=[])
            for k,v in self.photoz.items():
                photoz.attrs[k]=v


    def dumpPhoSim(self, fname, header="", manyFiles=False, sedName="../sky/sed_flat.txt", 
                   objtype="sersic2D", ssize=2.0*u.arcsec):
        """
        Writes out catalog in format that phosim can chew.

        Parameters
        ----------
        fname : string
                Filename to dump it to. If manyFiles is true, it must contain
                a formatting string for each file, e.g. name%04d.txt
        header : string, optional
                Write this at the begginng of every file.
        manyFiles: boolean, optional
                If true, create one file per object, otherwise create one file
                for the entire catalog.
        sedName : string to put where sedName goes into phosim file.
        objtype : string
                  at the moment just point and sersic2D are supported
        ssize  : astropy quantity
                 sersic size
        """

        if objtype not in ["point","sersic2D"]:
            print ("Bad obj type",objtype)
            stop()

        ssize=float(ssize/u.arcsec)

        if manyFiles:
            tosave=[self.data[i:i+1] for i in range(len(self.data))]
        else:
            tosave=[self.data]

        oc=0
        for cc,lines in enumerate(tosave):
            if manyFiles:
                of=open(fname%cc,'w')
            else:
                of=open(fname,'w')
            of.write(header)
            for obj in lines:
                of.write ("object {ID} {RA} {DEC} {MAG_NORM} {SED_NAME} "
                          "{REDSHIFT} {GAMMA1} {GAMMA2} {KAPPA} {DELTA_RA} "
                          "{DELTA_DEC} {OBJTYPE} ".format(ID=oc, RA=float(obj["ra"]*u.rad/u.deg),
                          DEC=float(obj["dec"]*u.rad/u.deg), MAG_NORM=obj["rmag"], SED_NAME=sedName,
                          REDSHIFT=obj["z"], GAMMA1=obj["g1"], GAMMA2=obj["g2"], KAPPA=0.0,
                                                DELTA_RA=0.0, DELTA_DEC=0.0, OBJTYPE=objtype))

                if objtype=="sersic2D":
                    s=Shear(g1=obj["e1"], g2=obj["e2"])
                    beta=s.beta.rad()/np.pi*180.
                    rat=np.exp(s.eta)

                    of.write ("{major} {minor} {beta} {sersic}".format(major=ssize/np.sqrt(rat), 
                                                    minor=ssize*np.sqrt(rat), beta=beta, sersic=1))
                
                of.write("\n")
                oc+=1
            of.close()
            
                            
        
