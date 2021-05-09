import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *
from nptdms import TdmsFile
from signalprocessing import *

plt.rcParams["font.size"] = 18


# Define Class TDMS
#  TDMS class has functions to load, initialize variables, view waterfall contour map

class Tdms:
    def __init__(self):
        self.figname=False
        self.filename=False
    

    def load_variables(self,filepath):
        self.file=TdmsFile(filepath)
        self.groupname=self.file.groups()[0]
        self.channels=self.file.group_channels(self.groupname)
        self.nchannels=len(self.channels)
        self.nt=len(self.channels[0].data)
        self.header=self.file.object()
        # self.file=TdmsFile(filepath)
        # self.groupname=self.file.groups()[0]
        # self.channels=self.groupname.channels()
        # self.nchannels=len(self.channels)
        # self.nt=len(self.channels[0].data)
        # self.header=self.file.object()        
        
        self.init_variables()

        
        self.data=np.ndarray((self.nt,self.nchannels))
        for i in range(self.nchannels):
            self.data[::,i]=self.channels[i].data
        self.info()

    def init_variables(self):
        self.dt=1.0/self.header.property("SamplingFrequency[Hz]")
        self.tini=0.0
        self.tend=(self.nt-1)*self.dt+self.tini
        self.tt=np.arange(self.tini,self.tend+self.dt,self.dt)
        self.tt=self.tt[0:self.nt]
        self.zini=self.header.property("Start Distance (m)")
        self.zend=self.header.property("Stop Distance (m)")
        self.dz=self.header.property("Fibre Length Multiplier")*self.header.property("SpatialResolution[m]")
        self.zz=np.arange(self.zini,self.zend,self.dz)
        self.zz=self.zz[0:self.nchannels]

    def info(self):              
        print("===================")
        print("Filename :",self.header.property('name'))
        print("nsamples :",self.nt)
        print("nchannels:",self.nchannels)
        print("Time(ini,end,int):",self.tini,self.tend,self.dt)
        print("Distance(ini,end,int):",self.zini,self.zend,self.dz)
        print("(max,min):",np.max(self.data),np.min(self.data))
        
    def print_headers(self):
        for name, value in self.header.properties.items():
            print("{0}: {1}".format(name,value))
            
    def view(self):
        clip=0.8
        vmin=np.min(self.data)*clip
        vmax=np.max(self.data)*clip
        
        nt=self.nt
        nchannel=self.nchannels
#        xx3,zz3=np.meshgrid(np.linspace(0,(nt),nt),np.linspace(0,(nchannel),nchannel,endpoint=True))
        xx3,zz3=np.meshgrid(self.tt,self.zz)
        fig=plt.figure(figsize=(6,4))
        ax=fig.add_subplot(1,1,1)
        cax=plt.pcolormesh(xx3,zz3,self.data[::,::].T,cmap='gray',vmin=vmin,vmax=vmax)
        ax.invert_yaxis()
        ax.set_xlabel('Time [sec]')
        ax.set_ylabel('Distance [m]')
        
    def view_sparse(self,subsampling=10):
            clip=0.8
            vmin=np.min(self.data)*clip
            vmax=np.max(self.data)*clip
            
            nt=self.nt
            nchannel=self.nchannels
    #        xx3,zz3=np.meshgrid(np.linspace(0,(nt),nt),np.linspace(0,(nchannel),nchannel,endpoint=True))
            xx3,zz3=np.meshgrid(self.tt[0:-1:subsampling],self.zz)
            fig=plt.figure(figsize=(6,4))
            ax=fig.add_subplot(1,1,1)
            cax=plt.pcolormesh(xx3,zz3,self.data[0:-1:subsampling,::].T,cmap='gray',vmin=vmin,vmax=vmax)
            ax.invert_yaxis()
            ax.set_xlabel('Time [sec]')
            ax.set_ylabel('Distance [m]')

def main():
    print("TDMS_Functions")
    # Set the data path to the TDMS files
    dirpath='D:\\Research\\OpticSensing\\pretest\\RITE_first_trial\\'


# In[227]:


# Example to use TDMS class and functions
#

# .. Declare the class Tdms as test1
    test1=Tdms()

# .. Load variables from the TDMS file
    test1.load_variables(dirpath+"test01_UTC_20190316_082052.162.tdms")


# In[228]:


# .. View the contour map
    test1.view()


# In[229]:


# .. View the all header variables
    test1.print_headers()
# In[226]:


if __name__== '__main__':
    main()

def load_data_with_filter(filepath, ar_denoise=(-90,-20), order=5, hicut=40, 
                          locut=2, depthinit=True, denoiseflg=True, bpfflg=True, 
                          sliceflg={'bc880':True, 'hwc250':True,'stc250':True,
                                    'stcSurf':True,'hwcSurf':True,'welma':False},
                          slicevalue={'bc880':(3938,4690),'hwc250':(1016,1299),
                                      'stc250':(3243,3494),'stcSurf':(2800,3243),
                                      'hwcSurf':(525,1016),'welma':(2916,3615)},
                          fileformat='tdms', temptdms=False, tdmschange=False,
                          nsamples=16000, nchannels=4992):
    
    def chunckrun(slx):
        if bpfflg:
            slx.data=apply_bpf(slx.data,locut,hicut,1000,order=5)
        if depthinit:
            init_zz(slx)
        
    event1023=Model()
    #filepath='/mnt/h/20200210/connected whole_UTC_20200212_103700.000.tdms'
    if fileformat=='tdms':
        event1023.entire=load_tdms_as_Fibers(filepath)
    elif fileformat=='segy':
        if tdmschange==True:
            event1023.entire=load_segy_as_Fibers_change(filepath, 
                                                        tdmsfile=temptdms,
                                                        nchannels=nchannels,
                                                        nsamples=nsamples)
        else:
            event1023.entire=load_segy_as_Fibers(filepath,tdmsfile=temptdms)
    
    print('Entire data size=',event1023.entire.data.shape)

    if denoiseflg:
        event1023.entire.data=denoising(event1023.entire,ar_denoise)
    
    if sliceflg['bc880']:
        event1023.bc880=slice_spacechunk(event1023.entire,slicevalue['bc880'])
        chunckrun(event1023.bc880)
    if sliceflg['stc250']:
        event1023.stc250=slice_spacechunk(event1023.entire,slicevalue['stc250'])
        chunckrun(event1023.stc250)
    if sliceflg['hwc250']:
        event1023.hwc250=slice_spacechunk(event1023.entire,slicevalue['hwc250'])
        chunckrun(event1023.hwc250)       
    if sliceflg['hwcSurf']:
        event1023.hwcSurf=slice_spacechunk(event1023.entire,slicevalue['hwcSurf'])
        chunckrun(event1023.hwcSurf)
    if sliceflg['stcSurf']:
        event1023.stcSurf=slice_spacechunk(event1023.entire,slicevalue['stcSurf'])
        chunckrun(event1023.stcSurf)
    if sliceflg['welma']:
        event1023.welma=slice_spacechunk(event1023.entire,slicevalue['welma'])
        chunckrun(event1023.welma)    
    return event1023

def load_segy_as_Fibers_change(segyfile,tdmsfile,nchannels=4992,nsamples=16000):
    slx2020=Tdms()
    slx2020.load_variables(tdmsfile)
    slx=Fibres()
    slx.nsamples=nsamples
    slx.nchannels=nchannels
    slx.nt=nsamples

    
    print('inside SEGY loading function',slx.nsamples,slx.nchannels)

    slx.zz=mycp(slx2020.zz[0:nchannels])
#    slx.data=mycp(slx2020.data)
    slx.read_data(segyfile,slx.nchannels,slx.nsamples,sgy=True,endian='big')
    print(segyfile)
    print('Inside, data shape=',slx.data.shape)
    slx.tt=mycp(slx2020.tt[0:nsamples])
    return slx    
