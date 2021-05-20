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
                  xlabel='Variable', ylabel='Channel', clabel='Data'):
  """
  Plot waterfall of anything

  INPUT:

  x: 1D array data at x-axis with shape (m,)
  y: 1D array data at y-axis with shape (n,). Usually y-axis is channels OR depths.
  Z: 2D array to be plotted with shape (n,m)

  OUTPUT:

  Waterfall plot of Z
  """
  plt.imshow(Z, aspect='auto', extent=(min(x), max(x), max(y), min(y)),
             cmap=cmap, vmin=vmin, vmax=vmax)
  h = plt.colorbar()
  h.set_label(clabel)
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
  h = plt.colorbar()
  h.set_label('Amplitude')
  plt.title('Spectrogram')
  plt.xlabel('Time [sec]')
  plt.ylabel('Frequency [Hz]')
  return Zxx.shape

def fk(data, fs_int=0.001, chan_int=1, plot=False, vels=None, cmap='jet', 
       vmin=None, vmax=None, xlim=(None,None), ylim=(None,None)):
  """
  Calculate FK (Frequency-Wavenumber) and plot

  INPUT:

  data: Data (2D array) has shape (m,n) where m: time samples, n: number of channels
  fs_int: Time samples interval (seconds). Default is 0.001
  chan_int: Channel interval (distance between 2 channels, m). Default is 1.
  plot: Option to plot FK. Default is False, no plot is created.
  vels: List of velocities to be plotted as contours

  OUTPUT:

  f: Frequencies (Hz)
  k: Wavenumbers
  sp: FK magnitude
  Plot of F-K domain
  """
  m, n = data.shape

  # Frequencies
  nt = m # Number of time samples
  nyq_f = nt//2
  f = np.fft.fftfreq(nt, d=fs_int)[slice(0,nyq_f)]
  
  # Wavenumbers
  nx = n # Number of channels
  nyq_k = nx//2
  k = np.fft.fftshift(np.fft.fftfreq(nx, d=chan_int))
  
  # Frequency-wavenumber power spectral density
  A = np.fft.fftshift(np.fft.fft2(data)[:,slice(0,nyq_f)],axes=0)
  sp2 = 2*(np.abs(A)**2) / (nt**2)
  sp2 = 10*np.log10(sp2)

  # Results
  k, sp2 = k[2:], sp2.T

  if plot==True:
    # Plot FK
    plt.imshow(sp2, extent=[max(k), min(k), min(f), max(f)], 
               aspect='auto',cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
    h = plt.colorbar()
    h.set_label('Power Spectra [dB]')

    if vels!=None:
      # Plot velocities contours
      for c in vels:
        plt.plot(k, k*c, '--', label='{} m/s'.format(c))
      plt.legend(fontsize=10)

    plt.title('F-K Plot')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Wavenumber [1/m]')
    plt.xlim(xlim)
    plt.ylim(ylim)
  return f, k, sp2
