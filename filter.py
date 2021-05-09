import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *
from nptdms import TdmsFile
from scipy.signal import butter, lfilter, filtfilt

def depth_stack(data,nstack=1):
    nsamples,nchannels=data.shape
    nchout=int(nchannels/nstack)
    dataout=np.ndarray((nsamples,nchout))
    for i in range(nchout):
        dataout[:,i]=np.sum(data[:,i:i+nstack],axis=1)/nstack
    return dataout

def apply_moveout(datain,moveout,dt=0.001):
    shift=-int(moveout[-1]/dt)
    dataout=np.zeros_like(datain)
#    ((datain.shape[0]-shift,len(moveout)))
    print(dataout.shape)
    for i in range(len(moveout)):
        dataout[shift::,i]=datain[shift+int(moveout[i]/dt):int(moveout[i]/dt)+datain.shape[0],i]
    
    return dataout

def get_meanvel(depth):
    velmodel=np.loadtxt('/mnt/d/Research/OpticSensing/ALJ.2017/ALJ.2017/velocity_vertical_Ichihara.txt')
    
    res1p=np.polyfit(velmodel[:,2], velmodel[:,3], 1)
    res2p=np.polyfit(velmodel[:,2], velmodel[:,3], 2)
    res3p=np.polyfit(velmodel[:,2], velmodel[:,3], 3)


    res3s=np.polyfit(velmodel[:,2], velmodel[:,4], 1)
    res3s=np.polyfit(velmodel[:,2], velmodel[:,4], 2)
    res3s=np.polyfit(velmodel[:,2], velmodel[:,4], 3)
    
    dep=np.arange(0,depth,1)
    vels=np.poly1d(res3s)(dep)
    velp=np.poly1d(res3p)(dep)
    
    return np.mean(velp),np.mean(vels)

def get_traveltime(meanvel,distance):
    return distance/meanvel


# def butter_bandpass(lowcut, highcut, fs, order=5):
#     nyq = 0.5 * fs
#     low = lowcut / nyq
#     high = highcut / nyq
#     b, a = butter(order, [low, high], btype='band')
#     return b, a


# def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
#     b, a = butter_bandpass(lowcut, highcut, fs, order=order)
#     y = lfilter(b, a, data)
#     return y

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def rms(data):
    nsample=data.shape[0]
    rms=np.sqrt(np.sum(data**2,axis=0)/nsample)
    return rms
def get_snr(data,time,ch,rmswindow=2):
    index1=np.argmax(data.tt>time)
    index2=np.argmax(data.tt>time-rmswindow)
    index3=np.argmax(data.tt>time+rmswindow)
    
    rmsnoise=rms(data.data[index2:index1,ch])
    rmssignal=rms(data.data[index1:index3,ch])
    
    snr=rmssignal/rmsnoise
    
    return snr
def get_snr_single(data,time,rmswindow=2):
    index1=np.argmax(data.tt>time)
    index2=np.argmax(data.tt>time-rmswindow)
    index3=np.argmax(data.tt>time+rmswindow)
    
    rmsnoise=rms(data.data[index2:index1])
    rmssignal=rms(data.data[index1:index3])
    
    snr=rmssignal/rmsnoise
    
    return snr
def rms_all(data):
    nsample=data.shape[0]
    rms=np.square(np.sum(data**2,axis=0)/nsample)
    return rms
def get_snr_all(data,time,rmswindow=2,depini=2834,depend=3664):
    index1=np.argmax(data.tt>time)
    index2=np.argmax(data.tt>time-rmswindow)
    index3=np.argmax(data.tt>time+rmswindow)
    
    rmsnoise=rms_all(data.data[index2:index1,depini:depend])
    rmssignal=rms_all(data.data[index1:index3,depini:depend])
    
    snr=rmssignal/rmsnoise
    
    return snr
    
def moving_average(data,nave):
    outdata=np.convolve(data,np.ones(nave)/float(nave),'same')
    return outdata
