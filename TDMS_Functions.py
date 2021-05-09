import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *
from nptdms import TdmsFile


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
