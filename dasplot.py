import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from signalprocessing import *
from scipy import signal

def plotWaveformTraces(event, magnify=1, lpf=None, color='black', alpha=0.2, 
                       figsize=(7,5), xlim=(None, None), ylim=(None, None), 
                       title="DAS Response"):
  """
  Plotting routine for DAS waveforms

  NOTE: X-axis is time, Y-axis is depth  

  INPUT:

  event: Event data (TDMS object)
  magnify: Scale of magnification. Default is 1, no magnification.
  lpf: Low-pass filter. 
    * Default is None: No filter is used
    * If used, specify as dictionary of f (low pass frequency) and fs (sampling
      frequency)
  xlim: Limit of time, tuple (start_time, end_time). Default is None, 
    data at all time samples will be plotted.          
  ylim: Limit of depth, tuple (start_depth, end_depth). Default is None, 
    data at all depths will be plotted.        
  """
  data = event.data
  t, z = event.tt, event.zz
  m, n = data.shape

  if ylim!=(None,None):
    # Limit of depth is passed
    id_loc = np.where((z>=ylim[0]) & (z<=ylim[1]))[0]
    id0, id1 = id_loc[0], id_loc[-1]
    data = data[:,id0:id1+1] # Tricky part

    z = z[id_loc]
    m, n = data.shape
    # return data

  plt.figure(figsize=figsize)
  for i, j in zip(range(n), reversed(range(n))):
    ampl = data[:,i]
    if lpf!=None:
      # Filter is used
      ampl = butter_lowpass_filter(ampl, lpf["f"], lpf["fs"], order=5)
    ampl = ampl + j * magnify
    if i==0:
      y0 = min(ampl) + ((max(ampl) - min(ampl)) / 2)
    if i==n-1:
      y1 = min(ampl) + ((max(ampl) - min(ampl)) / 2)

    plt.plot(t, ampl, color=color, alpha=alpha)  
  
  # Customize y axis tick
  tick0, tick1 = z[0], z[-1]
  y = np.linspace(y0, y1, 10)
  tick = np.linspace(tick0, tick1, 10)
  plt.yticks(y, np.round(tick))

  # Labels and limits
  plt.xlabel("Time [sec]")
  plt.ylabel("Depth [m]")
  plt.xlim(xlim)
  # plt.ylim(ylim)
  plt.title(title, pad=10)

def plotVSP(event, magnify=1, lpf=None, color='black', alpha=0.2, 
            figsize=(7,5), xlim=(None, None), ylim=(None, None), 
            title="DAS-VSP Response"):
  """
  Plotting DAS waveforms as VSP style

  NOTE: X-axis is depth, Y-axis is time

  INPUT:

  event: Event data (TDMS object)
  magnify: Scale of magnification. Default is 1, no magnification.
  lpf: Low-pass filter. 
    * Default is None: No filter is used
    * If used, specify as dictionary of f (low pass frequency) and fs (sampling
      frequency)
  xlim: Limit of depth, tuple (start_depth, end_depth). Default is None, 
    data at all depths will be plotted. 
  ylim: Limit of time, tuple (start_time, end_time). Default is None, 
    data at all time samples will be plotted.     
  """
  data = event.data
  t, z = event.tt, event.zz
  m, n = data.shape

  if xlim!=(None,None):
    # Limit of depth is passed
    id_loc = np.where((z>=xlim[0]) & (z<=xlim[1]))[0]
    id0, id1 = id_loc[0], id_loc[-1]
    data = data[:,id0:id1+1]

    z = z[id_loc]
    m, n = data.shape

  plt.figure(figsize=figsize)
  for i, j in zip(range(n), reversed(range(n))):
    ampl = data[:,i]
    if lpf!=None:
      # Filter is used
      ampl = butter_lowpass_filter(ampl, lpf["f"], lpf["fs"], order=5)
    ampl = (ampl + i * magnify)
    if i==0:
      x0 = min(ampl) + ((max(ampl) - min(ampl)) / 2)
    if i==n-1:
      x1 = min(ampl) + ((max(ampl) - min(ampl)) / 2)

    plt.plot(ampl, t, color=color, alpha=alpha)  
  
  # Customize y axis tick
  tick0, tick1 = z[0], z[-1]
  x = np.linspace(x0, x1, 10)
  tick = np.linspace(tick0, tick1, 10)
  plt.xticks(x, np.round(tick), rotation=45)

  # Labels and limits
  plt.ylim(ylim)
  plt.gca().invert_yaxis()
  plt.xlabel("Depth [m]")
  plt.ylabel("Time [sec]")
  plt.title(title, pad=10)

# Test plot
# lpf = None
# plotWaveformTraces(event.bc880, magnify=1e3, lpf=lpf, xlim=(13.8,15))

def plotWaterfall(x, y, Z, cmap='jet', vmin=None, vmax=None, 
                  xlim=(None,None), ylim=(None,None), title='Waterfall Plot',
                  xlabel='Variable', ylabel='Channel'):
  """
  Plot waterfall of anything

  INPUT:

  x: 1D array data at x-axis with shape (m,)
  y: 1D array data at y-axis with shape (n,). Usually y-axis is channels OR depths.
  Z: 2D array to be plotted with shape (n,m)

  OUTPUT:

  Waterfall plot of Z
  """
#   assert (len(x)==Z.shape[1]) & (len(y)==Z.shape[0]), 'Shape of data Z with x and y unmatched. Check each dimensions!'
  plt.imshow(Z, aspect='auto', extent=(min(x), max(x), max(y), min(y)),
             cmap=cmap, vmin=vmin, vmax=vmax)
  plt.colorbar()
  plt.xlim(xlim)
  plt.ylim(ylim)
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)

def stftSpectrogram(x, fs, cmap='jet', window='hann', nperseg=256, 
                    noverlap=None, nfft=None, detrend=False, 
                    return_onesided=True, boundary='zeros', padded=True):
  """
  Short Time Fourier Transform of a signal and plot Spectrogram

  INPUT:

  x: Trace data (1D array)
  fs: Sampling frequency (=1/sampling rate in sec)
  cmap: Matplotlib colormap

  window, nperseg, noverlap, nfft, detrend, return_onesided, boundary, and
  padded are inputs for "scipy.signal.stft" function. Read the docs for details:
  SciPy docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.stft.html

  OUTPUT:

  Plot of spectrogram
  """
  f, t, Zxx = signal.stft(x, fs=fs, window=window, nperseg=nperseg, 
                          noverlap=noverlap, nfft=nfft, detrend=detrend, 
                          return_onesided=return_onesided, boundary=boundary, 
                          padded=padded)
  plt.pcolormesh(t, f, np.abs(Zxx), cmap=cmap)
  plt.colorbar()
  plt.title('Spectrogram')
  plt.xlabel('Time [sec]')
  plt.ylabel('Frequency [Hz]')
  return Zxx.shape
