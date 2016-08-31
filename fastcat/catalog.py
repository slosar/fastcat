from __future__ import print_function, division

import window
import photoz
import numpy as np
import astropy.units as u
from numpy.lib import recfunctions

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
    
    See hdf5_format_doc for a more complete documentation

    On construction:
    Options
    -------
    N: int
       number of objects in the catalog
    meta: string
       string containing meta info
   """
    version=0.3
    
    def __init__ (self, N=0, fields=['ra','dec','z'],dNdz=None, bz=None,window=window.WindowBase(),
                  photoz=None,meta=None, read_from=None):
        if (read_from!=None):
            self.readH5(read_from)
        else:
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
        self.data=of["objects"].value
        self.meta=of["meta"].attrs
        if "dNdz" in of.keys():
            self.dNdz=of['dNdz'].value 
        if "bz" in of.keys():
            self.bz=of['bz'].value
        self.window=window.readWindowH5(of['window'])
        self.photoz=photoz.readPhotoZH5(of['photoz'])
        cversion=float(self.meta['version'])
        if cversion==0.1:
            print("updating 0.1 to version ", self.version)
            self.data=recfunctions.append_fields(self.data,'sigma_pz',(1+self.data["z_real_t"])*self.data["z_error"],
                                                 usemask=False)
            self.data=recfunctions.append_fields(self.data,'z',self.data["z_real_t"]+(1+self.data["z_real_t"])*self.data["z_error"],
                                                 usemask=False)
            self.data=self.data[ [ name for name in self.data.dtype.names if name not in ["z_real_t", "z_rsd_t","z_error"] ] ]
        if cversion==0.2:
            print("WARNING: upgrading from 0.2 to 0.3, photozs internally slightly inconsistent.")
            self.data=recfunctions.append_fields(self.data,'sigma_pz',(1+self.data["z"])*self.photoz.sigma,
                                                 usemask=False)

    def writeH5(self, fname):
        """ 
        Write Catalog to H5 file, specified as argument
        """
        of=h5py.File(fname, "w")
        if (self.meta):
            meta=of.create_dataset("meta",data=[])
            for v in self.meta.keys():
                meta.attrs[v]=self.meta[v]
            meta.attrs['version']=self.version
            
        dset=of.create_dataset("objects", data=self.data, chunks=True,
                               shuffle=True,compression="gzip", compression_opts=9)
        if type(self.dNdz)!=type(None):
            dset=of.create_dataset("dNdz", data=self.dNdz)
        if type(self.bz)!=type(None):
            dset=of.create_dataset("bz", data=self.bz)
        self.window.writeH5(of)
        pz=of.create_dataset("photoz",data=[])
        self.photoz.writeH5(pz)

    def setWindow(self,window,apply_to_data=True):
        """ 
        Sets window function to window
        and then optionally applies it to the current data by sampling probabilities
        """
        self.window=window
        if (apply_to_data):
            self.data=self.window.applyWindow(self.data)

    def setPhotoZ(self,photoz,apply_to_data=True):
        """ 
        Sets PZ description
        and then optionally applies it to the current data by sampling probabilities
        """
        self.photoz=photoz
        if (apply_to_data):
            self.data=self.photoz.applyPhotoZ(self.data)

    def appendCatalog(self, addcat):
        """
        Expands current catalog and then for those fields that are in common,
        add data from newcat.
        """
        N1=len(self.data)
        Nx=len(self.data.dtype)
        N2=len(addcat.data)
        newdata=np.zeros(((N1+N2),),dtype=self.data.dtype)
        newdata[0:N1]=self.data
        for n in addcat.data.dtype.names:
            if n in newdata.dtype.names:
                newdata[n][N1:]=addcat.data[n]
            else:
                print("Warning: not adding ",a," in catalog.appendCatalog")
        self.data=newdata


            
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
            
                            
        
