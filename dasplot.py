import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plotWaveformTraces(event, magnify=1, lpf=None, color='black', alpha=0.2, 
                       figsize=(7,5), xlim=(None, None), ylim=(None, None), 
                       title="DAS Response"):
  """
  Plotting routine for DAS waveforms

  INPUT:

  event: Event data (TDMS object)
  magnify: Scale of magnification. Default is 1, no magnification.
  lpf: Low-pass filter. 
    * Default is None: No filter is used
    * If used, specify as dictionary of f (low pass frequency) and fs (sampling
      frequency)
  """
  data = event.data
  t, z = event.tt, event.zz
  m, n = data.shape

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
  plt.ylabel("Distance [m]")
  plt.xlim(xlim)
  plt.ylim(ylim)
  plt.title(title, pad=10)

# Test plot
# lpf = None
# plotWaveformTraces(event.bc880, magnify=1e3, lpf=lpf, xlim=(13.8,15))
