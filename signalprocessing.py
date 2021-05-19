import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *
from nptdms import TdmsFile
from TDMS_Functions import *
from filter import *
import scipy.fft as ft
from scipy import interpolate, signal

def get_hhmmss(time):
    hh=int((time-np.mod(time,3600))/3600)
    time2=time-hh*3600
    mm=int((time2-np.mod(time2,60))/60)
    ss=time2-mm*60
    return hh,mm,ss

def get_seconds(time):
    hh,mm,ss=get_hhmmss_str(time)
    return seconds
def hhmmss_seconds(hh,mm,ss):
    seconds=hh*3600+mm*60+ss
    return seconds
def get_hours(time):
    hh,mm,ss=get_hhmmss_str(time)
    seconds=hhmmss_seconds(hh,mm,ss)
    hours=seconds/3600
    return hours
def get_hhmmss_str(time):
    hh,mm,ss=time.split(':')
    hh,mm,ss=(int(hh),int(mm),float(ss))
    return hh,mm,ss

def yymmdd(date,delimiter='/',flg='int'):
    tmp=date.split(delimiter)
    yy,mm,dd=(tmp[0],tmp[1],tmp[2])
    if flg=='int':
        yy,mm,dd=(int(yy),int(mm),int(dd))
    else:
        yy,mm,dd=yy.strip(),mm.strip(),dd.strip()
    return yy,mm,dd

def hhmmss(time,delimiter=':',flg='int'):
    if delimiter=='':
        hh=time[0:2]
        mm=time[2:4]
        ss=time[4:]
    else:
        tmp=time.split(delimiter)
        hh,mm,ss=(tmp[0],tmp[1],tmp[2])
    
    if flg=='int':
        hh,mm,ss=(int(hh),int(mm),int(float(ss)))
    elif flg=='float':
        hh,mm,ss=(float(hh),float(mm),float(ss))
    else:
        hh,mm,ss=hh.strip(),mm.strip(),ss.strip()
    return hh,mm,ss

def get_sdate_stime(filename):
    sdate,stime=filename[9:17],filename[18:28]
    return sdate,stime

def get_filesize(filename):
    return os.path.getsize(filename)
def get_nsamples(filesize,nchannels,flg='sgy',bytepersample=4):
    if flg=='sgy':
        localnsample=int(((filesize-3600)/nchannels-240)/bytepersample)
    else:
        localnsample=(filesize/nchannels/bytepersample)
    return localnsample
def add_time(ini,addtime):
    
    hh=ini[0]+addtime[0]
    mm=ini[1]+addtime[1]
    ss=ini[2]+addtime[2]
    
    if ss>=60:
        ss=ss-60
        mm=mm+1
    
    if mm>=60:
        mm=mm-60
        hh=hh+1
        
#    if hh>=24:
#        comment='Next day'
    
    return hh,mm,ss
    
class Model:
    def __init__(self):
#        print('Sccessfully defined')
        self.set_xlabel('Variables')
        self.figname=False
        self.filename=False
        self.header=np.zeros(120)        
        
    def set_filename(self,filename):
        self.filename=filename
        
    def set_xlabel(self,label):
        self.xlabel=label
    
    def set_data(self,data):
        self.data=copy.deepcopy(data)
        
        
    def calc_pressure(self):
        value=gravity*imega*self.data
        return value
        
    def calc_percentile(self):
        self.p10=np.percentile(self.data,10)
        self.p50=np.percentile(self.data,50)
        self.p90=np.percentile(self.data,90)
        
    def calc_randn(self, seed, mean, sigma,nrand):
        np.random.seed(seed)
        self.data=np.random.normal(mean,sigma,(nrand))
        
    def calc_invlognormal(self):
        data=(1.0/(1.0+np.exp(-1.0*self.data)))
        return data
        
class Tosoku:
    def __init__(self):
        self.figname=False
        self.filename=False
        self.EW=np.ndarray(0)
        self.NS=np.ndarray(0)
        self.VZ=np.ndarray(0)
    def load_variables(self,filepath):
        DUMMYFILE='/mnt/d/Research/OpticSensing/2019_03_Mobara/市原地表加速度計データ/'+'Dummy1.asc'
        num_headers=10
        self.header={'filepath':filepath}
        
        if os.path.exists(filepath):
            localfilepath=filepath
        else:
            localfilepath=DUMMYFILE
        
        f=open(localfilepath)
        for i in range(num_headers):
            line=f.readline().rstrip('\r\n')
            propertyname=line[:20]
            value=line[20:]
            self.header[propertyname]=value
        
        self.nsamples=int(self.header['Number of Data      '])
        self.nt=self.nsamples
        temp=self.header['Trigger Time        '].split()
        self.date=temp[0]
        self.time=temp[1]
        self.dt=1.0/float(self.header['Sampling Freq(Hz)   '])
        
        self.tini=0.0
        self.tend=(self.nt-1)*self.dt+self.tini
        self.tt=np.arange(self.tini,self.tend+self.dt,self.dt)
        
        tempdata=f.readlines()
        self.data=np.ndarray(self.nsamples)
        for i in range(self.nsamples):
            self.data[i]=tempdata[i].rstrip('\r\n')
            
        f.close
        
    def set_3comp(self,filehead):
        filepath1=filehead+'.01.asc'
        filepath2=filehead+'.02.asc'
        filepath3=filehead+'.03.asc'
        self.load_variables(filepath1);self.EW=self.data
        self.load_variables(filepath2);self.NS=self.data
        self.load_variables(filepath3);self.VZ=self.data
        
    def add_3comp(self,extra):
        self.EW=np.append(self.EW,extra.EW)
        self.NS=np.append(self.NS,extra.NS)
        self.VZ=np.append(self.VZ,extra.VZ)
    
    def process_day(self,day):
        for i in range(0,24):
            shour='{:02d}'.format(i)
            for j in range(0,6):
                smin='{:01d}'.format(j)
                localpathhead=datapath+day+'/2019'+day+shour+smin+'000.203.t3w'
#                print(localpathhead)
                tmp=Tosoku();tmp.set_3comp(localpathhead);self.add_3comp(tmp)
    def subsample(self):
        self.EWsub10=self.EW[::10]
        self.NSsub10=self.NS[::10]
        self.VZsub10=self.VZ[::10]
        self.EWsub100=self.EW[::100]
        self.NSsub100=self.NS[::100]
        self.VZsub100=self.VZ[::100]
    def export(self,filepath):
        self.EW.astype('float32').tofile(filepath+'EW')
        self.NS.astype('float32').tofile(filepath+'NS')
        self.VZ.astype('float32').tofile(filepath+'VZ')
    def reset_header_day(self):
        self.nt=len(self.EW)
        self.tini=0.0
        self.dt=1.0/100.0
        self.tend=(self.nt-1)*self.dt+self.tini
        self.tend10=(self.nt/10-1)*self.dt*10+self.tini
        self.tend100=(self.nt/100-1)*self.dt*100+self.tini
        self.tt=np.arange(self.tini,self.tend+self.dt,self.dt)
        self.ttsub10=np.arange(self.tini,self.tend10+self.dt*10,self.dt*10)
        self.ttsub100=np.arange(self.tini,self.tend100+self.dt*100,self.dt*100)
    
    def view_sub10(self,tlim,title="Title",lentime=60,absymin=0.1):
        tmin=tlim[0]
        tmax=tlim[1]
        plt.rcParams["font.size"] = 14
        ymin=-1*absymin;ymax=absymin
        index=np.where(self.VZsub10==0.0000)
        self.VZsub10[index]=np.nan
        self.NSsub10[index]=np.nan
        self.EWsub10[index]=np.nan
        
        self.tth=self.ttsub10/3600
        itmin=np.argmax(self.tth>=tmin)
        itmax=np.argmax(self.tth>=tmax)

        
        fig=plt.figure(figsize=(16,8))
        ax1=fig.add_subplot(3,1,1)
        ax1.set_title(title)
        ax1.plot(self.tth[itmin:itmax],zero_mean(self.EWsub10[itmin:itmax]))
#        ax.plot([time2,time2],[ymin,ymax],'r--')
        ax1.set_ylabel('EW')
        ax1.set_xlim(tmin,tmax)
        ax1.set_ylim(ymin,ymax)

        ax2=fig.add_subplot(3,1,2)
        ax2.plot(self.tth[itmin:itmax],zero_mean(self.NSsub10[itmin:itmax]))
#        ax.plot([time2,time2],[ymin,ymax],'r--')
        ax2.set_ylabel('NS')
        ax2.set_xlim(tmin,tmax)
        ax2.set_ylim(ymin,ymax)

        ax3=fig.add_subplot(3,1,3)
        ax3.plot(self.tth[itmin:itmax],zero_mean(self.VZsub10[itmin:itmax]))
#        ax3.plot([time2,time2],[ymin,ymax],'r--')
        ax3.set_ylabel('Vertical')
        ax3.set_xlabel('Hour')
        ax3.set_xlim(tmin,tmax)
        ax3.set_ylim(ymin,ymax)
        return (ax1,ax2,ax3,fig)
    def view(self,tlim,title="Title",lentime=60,absymin=0.1):
        tmin=tlim[0]
        tmax=tlim[1]
        plt.rcParams["font.size"] = 14
        ymin=-1*absymin;ymax=absymin
        index=np.where(self.VZ==0.0000)
        self.VZ[index]=np.nan
        self.NS[index]=np.nan
        self.EW[index]=np.nan
        
        self.tth=self.tt/3600
        itmin=np.argmax(self.tth>=tmin)
        itmax=np.argmax(self.tth>=tmax)

        
        fig=plt.figure(figsize=(16,8))
        ax1=fig.add_subplot(3,1,1)
        ax1.set_title(title)
        ax1.plot(self.tth[itmin:itmax],zero_mean(self.EW[itmin:itmax]))
#        ax.plot([time2,time2],[ymin,ymax],'r--')
        ax1.set_ylabel('EW')
        ax1.set_xlim(tmin,tmax)
        ax1.set_ylim(ymin,ymax)

        ax2=fig.add_subplot(3,1,2)
        ax2.plot(self.tth[itmin:itmax],zero_mean(self.NS[itmin:itmax]))
#        ax.plot([time2,time2],[ymin,ymax],'r--')
        ax2.set_ylabel('NS')
        ax2.set_xlim(tmin,tmax)
        ax2.set_ylim(ymin,ymax)

        ax3=fig.add_subplot(3,1,3)
        ax3.plot(self.tth[itmin:itmax],zero_mean(self.VZ[itmin:itmax]))
#        ax3.plot([time2,time2],[ymin,ymax],'r--')
        ax3.set_ylabel('Vertical')
        ax3.set_xlabel('Hour')
        ax3.set_xlim(tmin,tmax)
        ax3.set_ylim(ymin,ymax)
        return (ax1,ax2,ax3,fig)

    
    
def zero_mean(indata):
    meanval=np.nanmean(indata)
    outdata=indata-meanval
    return outdata

class Events():
    def __init__(self):
        self.file='false'
        
        
    def epidistance(self,event):
        R = 6373.0
        lat1 = radians(self.loc[0])
        lon1 = radians(self.loc[1])
        lat2 = radians(event.loc[0])
        lon2 = radians(event.loc[1])

        dlon = lon2 - lon1
        dlat = lat2 - lat1

        a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

        distance = R * c
        return distance
    def distance(self,event):
        epidist=self.epidistance(event)
        distance=np.sqrt(epidist**2+(self.dep-event.dep)**2)
        return distance
        
class Parosci:
    def __init__(self):
        self.figname=False
        self.filename=False
        self.EW=np.ndarray(0)
        self.NS=np.ndarray(0)
        self.VZ=np.ndarray(0)
        self.temp=np.ndarray(0)
        self.oneg=np.ndarray(0)
    def load_variables(self,filepath):
        DUMMYFILE='/mnt/d/Research/OpticSensing/2019_03_Mobara/市原地表加速度計データ/'+'Dummy1.asc'
        num_headers=0
        self.header={'filepath':filepath}
        
        if os.path.exists(filepath):
            localfilepath=filepath
        else:
            localfilepath=DUMMYFILE
        
        
        data = np.loadtxt(localfilepath,delimiter=",",
                         dtype=[('col1','S24'),('col2','S5'),('col3','f8'),('col4','f8'),
                                              ('col5','f8'),('col6','f8'),('col7','f8')])    
        nsample=data.size
        
        self.EW=np.ndarray(nsample)
        self.NS=np.ndarray(nsample)
        self.VZ=np.ndarray(nsample)
        self.temp=np.ndarray(nsample)
        self.oneg=np.ndarray(nsample)
        
        for i in range(nsample):
            self.EW[i]=data[i][2]
            self.NS[i]=data[i][3]
            self.VZ[i]=data[i][4]
            self.temp[i]=data[i][5]
            self.oneg[i]=data[i][6]
        
        self.nsamples=nsample
        self.nt=self.nsamples
        self.dt=0.1
        
        self.tini=0.0
        self.tend=(self.nt-1)*self.dt+self.tini
        self.tt=np.arange(self.tini,self.tend+self.dt,self.dt)
        
    def add_3comp(self,extra):
        self.EW=np.append(self.EW,extra.EW)
        self.NS=np.append(self.NS,extra.NS)
        self.VZ=np.append(self.VZ,extra.VZ)
    
    def add_temp1G(self,extra):
        self.temp=np.append(self.temp,extra.temp)
        self.oneg=np.append(self.oneg,extra.oneg)
    
    def process_day(self,day='03-26'):
        for i in range(0,24):
            shour='{:02d}'.format(i)
#            for j in range(0,6):
#                smin='{:01d}'.format(j)
            localpathhead=datapath+'2019-'+day+'_'+shour+'_log.txt'
#                print(localpathhead)
            tmp=Parosci();tmp.load_variables(localpathhead);self.add_3comp(tmp),self.add_temp1G(tmp)
    def subsample(self):
        self.EWsub10=self.EW[::10]
        self.NSsub10=self.NS[::10]
        self.VZsub10=self.VZ[::10]
        self.tempsub10=self.temp[::10]
        self.onegsub10=self.oneg[::10]
        self.EWsub100=self.EW[::100]
        self.NSsub100=self.NS[::100]
        self.VZsub100=self.VZ[::100]
        self.tempsub100=self.temp[::100]
        self.onegsub100=self.oneg[::100]
        self.EWsub500=self.EW[::500]
        self.NSsub500=self.NS[::500]
        self.VZsub500=self.VZ[::500]
        self.tempsub500=self.temp[::500]
        self.onegsub500=self.oneg[::500]
        self.EWsub1000=self.EW[::1000]
        self.NSsub1000=self.NS[::1000]
        self.VZsub1000=self.VZ[::1000]
        self.tempsub1000=self.temp[::1000]
        self.onegsub1000=self.oneg[::1000]
    def export(self,filepath):
        self.EW.astype('float32').tofile(filepath+'EW')
        self.NS.astype('float32').tofile(filepath+'NS')
        self.VZ.astype('float32').tofile(filepath+'VZ')
        self.temp.astype('float32').tofile(filepath+'TEMP')
        self.oneg.astype('float32').tofile(filepath+'1G')
    def reset_header_day(self):
        self.nt=len(self.EW)
        self.tini=0.0
        self.dt=1.0/10.0
        self.tend=(self.nt-1)*self.dt+self.tini
        self.tend10=(self.nt/10-1)*self.dt*10+self.tini
        self.tend100=(self.nt/100-1)*self.dt*100+self.tini
        self.tend500=(self.nt/500-1)*self.dt*500+self.tini
        self.tend1000=(self.nt/1000-1)*self.dt*1000+self.tini
        self.tt=np.arange(self.tini,self.tend+self.dt,self.dt)
        self.ttsub10=np.arange(self.tini,self.tend10+self.dt*10,self.dt*10)
        self.ttsub100=np.arange(self.tini,self.tend100+self.dt*100,self.dt*100)
        self.ttsub500=np.arange(self.tini,self.tend500+self.dt*500,self.dt*500)
        self.ttsub1000=np.arange(self.tini,self.tend1000+self.dt*1000,self.dt*1000)
    def view_sub10(self,tlim,title="Title",lentime=60,absymin=0.1):
        tmin=tlim[0]
        tmax=tlim[1]
        plt.rcParams["font.size"] = 14
        ymin=-1*absymin;ymax=absymin
        index=np.where(self.VZsub10==0.0000)
        self.VZsub10[index]=np.nan
        self.NSsub10[index]=np.nan
        self.EWsub10[index]=np.nan
        
        self.tth=self.ttsub10/3600
        itmin=np.argmax(self.tth>=tmin)
        itmax=np.argmax(self.tth>=tmax)

        
        fig=plt.figure(figsize=(16,8))
        ax1=fig.add_subplot(3,1,1)
        ax1.set_title(title)
        ax1.plot(self.tth[itmin:itmax],zero_mean(self.EWsub10[itmin:itmax]))
#        ax.plot([time2,time2],[ymin,ymax],'r--')
        ax1.set_ylabel('EW')
        ax1.set_xlim(tmin,tmax)
        ax1.set_ylim(ymin,ymax)

        ax2=fig.add_subplot(3,1,2)
        ax2.plot(self.tth[itmin:itmax],zero_mean(self.NSsub10[itmin:itmax]))
#        ax.plot([time2,time2],[ymin,ymax],'r--')
        ax2.set_ylabel('NS')
        ax2.set_xlim(tmin,tmax)
        ax2.set_ylim(ymin,ymax)

        ax3=fig.add_subplot(3,1,3)
        ax3.plot(self.tth[itmin:itmax],zero_mean(self.VZsub10[itmin:itmax]))
        ax3.plot([time2,time2],[ymin,ymax],'r--')
        ax3.set_ylabel('Vertical')
        ax3.set_xlabel('Hour')
        ax3.set_xlim(tmin,tmax)
        ax3.set_ylim(ymin,ymax)
        return (ax1,ax2,ax3,fig)
    def view(self,tlim,title="Title",lentime=60,absymin=0.1):
        tmin=tlim[0]
        tmax=tlim[1]
        plt.rcParams["font.size"] = 14
        ymin=-1*absymin;ymax=absymin
        index=np.where(self.VZ==0.0000)
        self.VZ[index]=np.nan
        self.NS[index]=np.nan
        self.EW[index]=np.nan
        
        self.tth=self.tt/3600
        itmin=np.argmax(self.tth>=tmin)
        itmax=np.argmax(self.tth>=tmax)

        
        fig=plt.figure(figsize=(16,8))
        ax1=fig.add_subplot(3,1,1)
        ax1.set_title(title)
        ax1.plot(self.tth[itmin:itmax],zero_mean(self.EW[itmin:itmax]))
#        ax.plot([time2,time2],[ymin,ymax],'r--')
        ax1.set_ylabel('EW')
        ax1.set_xlim(tmin,tmax)
        ax1.set_ylim(ymin,ymax)

        ax2=fig.add_subplot(3,1,2)
        ax2.plot(self.tth[itmin:itmax],zero_mean(self.NS[itmin:itmax]))
#        ax.plot([time2,time2],[ymin,ymax],'r--')
        ax2.set_ylabel('NS')
        ax2.set_xlim(tmin,tmax)
        ax2.set_ylim(ymin,ymax)

        ax3=fig.add_subplot(3,1,3)
        ax3.plot(self.tth[itmin:itmax],zero_mean(self.VZ[itmin:itmax]))
#        ax3.plot([time2,time2],[ymin,ymax],'r--')
        ax3.set_ylabel('Vertical')
        ax3.set_xlabel('Hour')
        ax3.set_xlim(tmin,tmax)
        ax3.set_ylim(ymin,ymax)
        return (ax1,ax2,ax3,fig)

    
segystruct1='iiiiiiihhh'+'hiiiiiiiih'+'hiiiihhhhh'+'hhhhhhhhhh' # 
segystruct2='hhhhhhhhhh'+'hhhhhhhhhh'+'hhhhhhhhhh'+'hiiiiihhii'
segystruct3='ihhhhhhhhi'+'i'

segystruct=segystruct1+segystruct2+segystruct3

class Fibres:
    def __init__(self):
#        print('Sccessfully defined')
        self.figname=False
        self.filename=False

    def load_data(self,filename,nchannel,nsample):
        test=Model()
        test.set_filename(filename)
        fd=open(test.filename,'rb')
        data=np.fromfile(fd,dtype = np.dtype('float32'))
        self.data=np.reshape(data,(nchannel,nsample))
        self.data=self.data.T# Outer mo
        self.nsamples=nsample
        self.nchannels=nchannel
    def load_data_tdms(self,filepath):
        self.file=TdmsFile(filepath)
        self.groupname=self.file.groups()[0]
        self.channels=self.file.group_channels(self.groupname)
        self.nchannels=len(self.channels)
        self.nt=len(self.channels[0].data)
        self.header=self.file.object()
        
#        self.init_variables()

        
        self.data=np.ndarray((self.nt,self.nchannels))
        for i in range(self.nchannels):
            self.data[::,i]=self.channels[i].data
#        self.info()
    def read_data(self,filename,nchannel=15700,nsample=30000,sgy=False,endian='big'):
        bytestrh=240
        bytesfloats=4
        bytesblock=nsample*bytesfloats+bytestrh # =6244
        fin=open(filename,'rb')
        if sgy==True:
            fin.seek(3600)
        tmpval=[Model() for x in range(nchannel)]
        for i in range(nchannel):
            tmp=fin.read(bytestrh)
            if endian=='big':
                tmpval[i]=np.fromfile(fin,dtype='>f',count=nsample)
            else:
                tmpval[i]=np.fromfile(fin,dtype='<f',count=nsample)
        fin.close()
        self.data=np.array(tmpval).T
        self.nsamples=nsample
        self.nchannels=nchannel
        #return val

    def init_para(self,temp):
        self.dt=temp.dt
        self.tini=temp.tini
        self.tend=self.tini+(self.nsamples-1)*self.dt
        self.tt=np.arange(self.tini,self.tend+self.dt,self.dt)
        self.tt=self.tt[0:self.nsamples]
        self.nt=self.nsamples
        self.zini=temp.zini
        self.dz=temp.dz
        self.zend=self.zini+(self.nchannels-1)*self.dz
        self.zz=np.arange(self.zini,self.zend+self.dz,self.dz)
        self.zz=self.zz[0:self.nchannels]
    def view_sparse_as(self,ax,subsampling=10,clip=0.01,color='gray',myabs=-1):
#        clip=0.01
        if myabs<0:
            scaler=np.maximum(np.abs(np.min(self.data)),np.abs(np.max(self.data)))
            vmin=-scaler*clip
            vmax=scaler*clip
        else:
            vmin=-myabs
            vmax=myabs
        
        nt=self.nt
        nchannel=self.nchannels
#        xx3,zz3=np.meshgrid(np.linspace(0,(nt),nt),np.linspace(0,(nchannel),nchannel,endpoint=True))
        xx3,zz3=np.meshgrid(self.tt[0:-1:subsampling],self.zz)
#        fig=plt.figure(figsize=(6,4))
#        ax=fig.add_subplot(1,1,1)
        cax=plt.pcolormesh(xx3,zz3,self.data[0:-1:subsampling,::].T,cmap=color,vmin=vmin,vmax=vmax)
        ax.invert_yaxis()
        ax.set_xlabel('Time [sec]')
        ax.set_ylabel('Distance [m]')
        return ax
    def view_sparse_as_vertical(self,ax,subsampling=10,clip=0.01,color='gray',myabs=-1):
#        clip=0.01
        if myabs<0:
            scaler=np.maximum(np.abs(np.min(self.data)),np.abs(np.max(self.data)))
            vmin=-scaler*clip
            vmax=scaler*clip
        else:
            vmin=-myabs
            vmax=myabs
        
        nt=self.nt
        nchannel=self.nchannels
#        xx3,zz3=np.meshgrid(np.linspace(0,(nt),nt),np.linspace(0,(nchannel),nchannel,endpoint=True))
        zz3,xx3=np.meshgrid(self.zz,self.tt[0:-1:subsampling])
#        fig=plt.figure(figsize=(6,4))
#        ax=fig.add_subplot(1,1,1)
        cax=plt.pcolormesh(zz3,xx3,self.data[0:-1:subsampling,::],cmap=color,vmin=vmin,vmax=vmax)
        ax.invert_yaxis()
        ax.set_ylabel('Time [sec]')
        ax.set_xlabel('Distance [m]')
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        return ax
    def set_zz(self,base,index):
        self.zz=np.linspace(base.zz[index],base.zz[index+self.nchannels],self.nchannels,endpoint=True)
    def set_para(self,dt,tini,dz,zini):
        self.dt=dt
        self.tini=tini
        self.tend=self.tini+(self.nsamples-1)*self.dt
        self.tt=np.arange(self.tini,self.tend+self.dt,self.dt)
        self.tt=self.tt[0:self.nsamples]
        self.nt=self.nsamples
        self.zini=zini
        self.dz=dz
        self.zend=self.zini+(self.nchannels-1)*self.dz
        self.zz=np.arange(self.zini,self.zend+self.dz,self.dz)
        self.zz=self.zz[0:self.nchannels]        


def shift_utc(date,time,shift=-9):
    year,month,day=date[0],date[1],date[2]
    hh,mm,ss=time[0],time[1],time[2]
    nhh=hh+shift
        
    if nhh<0:
        nday=day-1
        nhh=nhh+24
    else:
        nday=day
        
    if nday<0:
        nmonth=month-1
        if nmonth in [1,3,5,7,8,10,12]:
            nday=nday+31
        elif nmonth==2:
            nday=nday+28
        else:
            nday=nday+30
    else:
        nmonth=month
        
    return (year,nmonth,nday),(nhh,mm,ss)

def shift_time(date,time,shift):
    year,month,day=date[0],date[1],date[2]
    hh,mm,ss=time[0],time[1],time[2]
    sfhh,sfmm,sfss=shift[0],shift[1],shift[2]
    
    nhh,nmm,nss=hh+sfhh,mm+sfmm,ss+sfss
        
    if nss<0:
        nmm=nmm-1
        nss=nss+60
    elif nss>60:
        nmm=nmm+1
        nss=nss-60
        
    if nmm<0:
        nhh=nhh-1
        nmm=nmm+60
    elif nmm>60:
        nhh=nhh+1
        nmm=nmm-60
        
    if nhh<0:
        nday=day-1
        nhh=nhh+24
    elif nhh>24:
        nhh=nhh-24
        nday=day+1
    else:
        nday=day
        
    if nday<0:
        nmonth=month-1
        if nmonth in [1,3,5,7,8,10,12]:
            nday=nday+31
        elif nmonth==2:
            nday=nday+28
        else:
            nday=nday+30
    else:
        nmonth=month
        
    return (year,nmonth,nday),(nhh,nmm,nss)

def copy_Fibres(old):
    import copy
    slx=copy.deepcopy(old)
    return slx
def slice_spacechunk(slxin,ar_dep):
    ch1=np.argmax(slxin.zz>=ar_dep[0])
    ch2=np.argmax(slxin.zz>=ar_dep[1])
    slx=copy_Fibres(slxin)
    slx.nchannels=ch2-ch1
    slx.zz=np.zeros(slx.nchannels)
    slx.zz=slxin.zz[ch1:ch2]
    slx.data=slxin.data[::,ch1:ch2]
    return slx
def slice_timechunk(slxin,ar_tt):
    ch1=np.argmax(slxin.tt>=ar_tt[0])
    ch2=np.argmax(slxin.tt>=ar_tt[1])
    slx=copy_Fibres(slxin)
    slx.nsamples=int(ch2-ch1)
    slx.tt=np.zeros(slx.nsamples)
    slx.tt=slxin.tt[ch1:ch2]
    slx.data=slxin.data[ch1:ch2,::]
    return slx

def denoising(impactor2020,ar_dep):
    index1=np.argmax(impactor2020.zz>ar_dep[0])
    index2=np.argmax(impactor2020.zz>ar_dep[1])
    noise=np.average(impactor2020.data[::,index1:index2],axis=1)
    stripenoise=np.zeros_like(impactor2020.data)
    size1,size2=impactor2020.data.shape
    for i in range(size2):
        stripenoise[::,i]=noise[:]
    return impactor2020.data-stripenoise
def denoising_out(impactor2020,noise,ar_dep):
    index1=np.argmax(noise.zz>ar_dep[0])
    index2=np.argmax(noise.zz>ar_dep[1])
    noise=np.average(noise.data[::,index1:index2],axis=1)
    stripenoise=np.zeros_like(impactor2020.data)
    size1,size2=impactor2020.data.shape
    for i in range(size2):
        stripenoise[::,i]=noise[:]
    return impactor2020.data-stripenoise




def mycp(obj):
    import copy
    return copy.deepcopy(obj)

# def load_tdms_as_Fibers(outpath):
#     slx2020=Tdms()
#     slx2020.load_variables(outpath)
#     slx=Fibres()
#     slx.nsamples=mycp(slx2020.nt)
#     slx.nchannels=mycp(slx2020.nchannels)
#     slx.nt=mycp(slx2020.nt)


#     slx.zz=mycp(slx2020.zz)
#     slx.data=mycp(slx2020.data)
#     slx.tt=mycp(slx2020.tt)
#     return slx

def load_segy_as_Fibers(segyfile,tdmsfile):
    slx2020=Tdms()
    slx2020.load_variables(tdmsfile)
    slx=Fibres()
    slx.nsamples=mycp(slx2020.nt)
    slx.nchannels=mycp(slx2020.nchannels)
    slx.nt=mycp(slx2020.nt)

    
    print('inside SEGY loading function',slx.nsamples,slx.nchannels)

    slx.zz=mycp(slx2020.zz)
#    slx.data=mycp(slx2020.data)
    slx.read_data(segyfile,slx.nchannels,slx.nsamples,sgy=True,endian='big')
    print(segyfile)
    print('Inside, data shape=',slx.data.shape)
    slx.tt=mycp(slx2020.tt)
    return slx

def init_zz(slx):
    slx.zz=slx.zz-slx.zz[0]

def apply_bpf(indata,lowcut,highcut,fs,order):
    len1,len2=indata.shape
    outdata=np.zeros_like(indata)
    for i in range(len2):
        outdata[:,i]=butter_bandpass_filter(indata[:,i],lowcut,highcut,fs,order)
        
    return outdata
    
def apply_depth_correction(slx,helicalcorr=0.881505039883044):
    slx.zz=slx.zz*helicalcorr  # sws3ST[1].zz[-1]/sws3ST[1].zz[-1]
    
def bc880_event(filepath,ar_denoise=(-90,-20),order=5,hicut=40,locut=2,
                depthinit=True,denoiseflg=True,bpfflg=True,
                fileformat='tdms',temptdms=False):
    event1023=Model()
    #filepath='/mnt/h/20200210/connected whole_UTC_20200212_103700.000.tdms'
    if fileformat=='tdms':
        event1023.entire=load_tdms_as_Fibers(filepath)
        print(event1023.entire.data.shape)
    elif fileformat=='segy':
        event1023.entire=load_segy_as_Fibers(filepath,tdmsfile=temptdms)
    if denoiseflg:
        event1023.entire.data=denoising(event1023.entire,ar_denoise)
    
#    event1023.hwc500=slice_spacechunk(event1023.entire,hwc500)
#    event1023.stc500=slice_spacechunk(event1023.entire,stc500)
#    event1023.hwc250=slice_spacechunk(event1023.entire,hwc250)
#    event1023.stc250=slice_spacechunk(event1023.entire,stc250)
    event1023.bc880=slice_spacechunk(event1023.entire,bc880)
#    event1023.hwcSurf=slice_spacechunk(event1023.entire,hwcSurf)
#    event1023.stcSurf=slice_spacechunk(event1023.entire,stcSurf)
    if bpfflg:
        event1023.bc880.data=apply_bpf(event1023.bc880.data,locut,hicut,1000,order=5)
#    event1023.hwc500.data=apply_bpf(event1023.hwc500.data,locut,hicut,1000,order=5)
#    event1023.hwc250.data=apply_bpf(event1023.hwc250.data,locut,hicut,1000,order=5)
#    event1023.stc500.data=apply_bpf(event1023.stc500.data,locut,hicut,1000,order=5)
#    event1023.stc250.data=apply_bpf(event1023.stc250.data,locut,hicut,1000,order=5)
#    event1023.hwcSurf.data=apply_bpf(event1023.hwcSurf.data,locut,hicut,1000,order=5)
#    event1023.stcSurf.data=apply_bpf(event1023.stcSurf.data,locut,hicut,1000,order=5)
    
    if depthinit:
        init_zz(event1023.bc880)
#        init_zz(event1023.hwc250)
#        apply_depth_correction(event1023.hwc250)
#        init_zz(event1023.stc250)
#    init_zz(event1023.bc880)
    
    return event1023

def normalize(x, method='MinMaxScaler'):
  """
  Normalize waveform

  INPUT:

  x: Data (2D array)
  method: Method options: 'MinMaxScaler', 'LogNorm', 'PowerTransformer', 'StandardScaler'

  OUTPUT:

  Normalized 2D array
  """
  if method=='MinMaxScaler':
    # Note that most min-max scaling algorithms (including Sklearn) only
    # scale from 0 to 1. The following has been adapted for negative input 
    # such that it would scale from -1 to 1. 
    xmax = np.amax(np.abs(x))
    range = 1
    return np.array([[float(val) / xmax * range for val in row] for row in x])
  if method=='LogNorm':
    return np.array([[np.log10(val) if val>0 else -np.log10(np.abs(val)) for val in row] for row in x])
  else:
    # Scikit-Learn transformers
    from sklearn.preprocessing import PowerTransformer, StandardScaler
    if method=='StandardScaler':
      scale = StandardScaler()
    if method=='PowerTransformer':
      scale = PowerTransformer()
    return scale.fit_transform(x)


def difference(waveform1, waveform2):
  """
  Calculate difference between two waveforms

  INPUT:

  waveform: Waveform data. Can be used to compare two different types of 
    DAS, here the STC (straight cable), HWC (helically wound cable), 
    or BC (behind casing).
    NOTE: waveform 1 must be inside waveform 2. 
  
  OUTPUT:

  dif: Waveform difference
  """
  data1, data2 = waveform1.data, waveform2.data
  x, y = waveform2.zz, waveform2.tt
  xnew, ynew = waveform1.zz, waveform1.tt
  f = interpolate.interp2d(x, y, data2)

  # Interpolate event 2 using (t,z) of event 1
  data2_i = f(xnew, ynew)

  # Now both events have same shape, take difference
  dif = data2_i - data1
  return dif

def fftSpectrum(x, fs, window=1, plot=True, flim=(None,None)): 
  """
  Calculate and amplitude spectrum
  
  INPUT:
  
  x: Trace data (1D array)
  fs: Sampling rate (milliseconds)
  window: Rolling window to smooth spectrum. Default is 1 (no smoothening)
  plot: Option to plot spectrum. Default is True.
  flim: Range of frequencies to plot. Default is None (up to Nyquist frequency)
  
  OUTPUT:
  
  frqs: Frequencies (1D array, Hz)
  frqAmp: Amplitude (1D array)
  Plot of spectrum
  """
  spectrum = ft.fft(x)
  frqBins = int(spectrum.size/2)
  frqAmp = np.absolute(spectrum[:frqBins]) 

  # Rolling average to smooth spectrum
  rolling = lambda x, w: np.convolve(x, np.ones(w), 'same') / w
  frqAmp = rolling(frqAmp, window)
  
  # Frequencies of interest
  NyquistFrq = fs/2.0 # the Nyquist frequency
  frqs = np.linspace(0, NyquistFrq, num=frqBins)
  
  if plot==True:
#       plt.figure()
      plt.plot(frqs, frqAmp, 'r')
      plt.xlabel('Frequency [Hz]')
      plt.ylabel('Amplitude')
      plt.xlim(flim)
      plt.title('Amplitude Spectrum')
  return frqs, frqAmp

# def stftSpectrogram(x, fs, cmap='jet', window='hann', nperseg=256, 
#                     noverlap=None, nfft=None, detrend=False, 
#                     return_onesided=True, boundary='zeros', padded=True):
#   """
#   Short Time Fourier Transform of a signal and plot Spectrogram

#   INPUT:

#   x: Trace data (1D array)
#   fs: Sampling frequency (=1/sampling rate in sec)
#   cmap: Matplotlib colormap

#   window, nperseg, noverlap, nfft, detrend, return_onesided, boundary, and
#   padded are inputs for "scipy.signal.stft" function. Read the docs for details:
#   SciPy docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.stft.html

#   OUTPUT:

#   Plot of spectrogram
#   """
#   f, t, Zxx = signal.stft(x, fs=fs, window=window, nperseg=nperseg, 
#                           noverlap=noverlap, nfft=nfft, detrend=detrend, 
#                           return_onesided=return_onesided, boundary=boundary, 
#                           padded=padded)
#   plt.pcolormesh(t, f, np.abs(Zxx), cmap=cmap)
#   plt.colorbar()
#   plt.title('Spectrogram')
#   plt.xlabel('Time [sec]')
#   plt.ylabel('Frequency [Hz]')
#   return Zxx.shape
