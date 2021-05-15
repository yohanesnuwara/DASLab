import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
