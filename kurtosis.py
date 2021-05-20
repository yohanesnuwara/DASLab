import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from filter import *

def kurtosis(data, window=100):
  """ 
  Calculate kurtosis over a specified window

  INPUT:

  data: Data (in seismic context: amplitudes), 1D array
  window: Kurtosis calculation window. Default is 100.

  OUTPUT:

  kurt: Kurtosis value, 1D array
  """
  kurt = []
  for i in range(len(data)):
    a = data[i:i+window]
    std, mean = np.std(data), np.mean(data)
    y = np.sum((a - mean)**4) / window
    k = y / std**4
    kurt.append(k)
  return kurt

def cutTrace(event, no_trace, cut_trace):
  """
  Cut trace given a specified window

  INPUT:

  event: Event data (TDMS object)
  no_trace: Trace number
  cut_trace: Cut trace window. Specify as tuple (start,end) with start is the
    start time and end is the end time

  OUTPUT:

  tcut: Trace time samples after cut operation
  ycut: Amplitude after cut operation
  """
  # Cut trace 
  t, z = event.tt, event.zz # axes
  index_cut = [int(np.where(t==cut_trace[0])[0]), int(np.where(t==cut_trace[1])[0])]
  tcut = t[index_cut[0]:index_cut[1]]

  fibre_data = event.data
  y = fibre_data[:,no_trace] 
  ycut = y[index_cut[0]:index_cut[1]]  
  return tcut, ycut

def kurtosisFindArrival(event, no_trace, cut_trace, win, lpf=None,
                          window=100, plot=False, pick_color=None, title=None):
  """
  Autopicking arrival times using Kurtosis and plot

  INPUT:

  event: Event data (TDMS object)
  cut_trace: Cut trace window. Specify as tuple (start,end) with start is the
    start time and end is the end time
  win: Searching window for time arrivals. It is a LIST of TUPLES, for instance
    [(x0,y0), (x1,y1), ..., (xn, yn)] where n is the number of arrivals you want
    to find, xn is the the starting range, and yn is the ending range. In case
    you want to find P and S arrival, then your win will be [(xp,yp), (xs, ys)]
  lpf: Low-pass filter. 
    * Default is None: No filter is used
    * If used, specify as dictionary of f (low pass frequency) and fs (sampling
      frequency)
  window: Kurtosis calculation window. Default is 100.
  plot: Option to plot the result. Default is False, so no plot is produced, but
    the following OUTPUTS will be produced.

  OUTPUT:

  t: Autopicked arrival times 
  A: Autopicked arrival kurtosis

  NOTE: Both t and A are LIST, with size that depends on your specified "win"
        For two arrivals i.e. P and S arrivals, LIST will have shape (2,)
  """  
  # Cut trace 
  tcut, ycut = cutTrace(event, no_trace, cut_trace)

  if lpf!=None:
    # Low pass filter is applied
    f, fs = lpf["f"], lpf["fs"]
    ycut = butter_lowpass_filter(ycut, f, fs, order=5)

  # Kurtosis calculation
  ycut_kurt = kurtosis(ycut, window=window)  
  
  assert type(win) is list, "Your window must be a list. Example: [(10,20)] if only one window, or [(10,20),(60,70)] if two windows."
  
  t, A = [], []
  for i in range(len(win)):
    # For every window
    win1 = np.where((tcut >= win[i][0]) & (tcut <= win[i][1]))[0]
    win1_start, win1_end = win1[0], win1[-1]

    win1_trace = ycut_kurt[win1_start:win1_end]

    # Search the highest point
    max_win1 = max(win1_trace)
    Ap = max_win1

    # Search time at the highest point (arrivals)
    tp = np.where(ycut_kurt == max_win1)[0]
    tp = tcut[tp][0]

    t.append(tp)
    A.append(Ap)

  if plot==True:
    # Plot seismogram and kurtosis
    plt.figure(figsize=(11,8))

    plt.subplot(2,1,1)
    plt.plot(tcut, ycut, color="black")
    plt.title("Trace Cut No. {} of {}".format(no_trace, title))
    plt.xlabel("Time [sec]")
    plt.ylabel("Amplitude")
    plt.xlim(min(tcut), max(tcut))
    plt.grid()

    # Plot autopicked arrivals
    for i in range(len(t)):
      if pick_color==None:
        plt.axvline(t[i])
      if pick_color!=None:
        plt.axvline(t[i], color=pick_color[i])

    plt.subplot(2,1,2)
    plt.title("Kurtosis of Trace No. {}".format(no_trace))
    plt.plot(tcut, ycut_kurt, color="red")
    plt.xlabel("Time [sec]")
    plt.ylabel("Kurtosis")
    plt.xlim(min(tcut), max(tcut))
    plt.grid()

    plt.tight_layout(1.1)
    plt.show()

  return t, A

# def kurtosisHeatmap(kurt, time_axis, n_traces, cmap='jet', figsize=(7,5),
#                     vmin=None, vmax=None):
#   """
#   Plot heatmap of kurtosis of all traces

#   INPUT:

#   kurt: 2D array of calculated kurtosis (output of running "pickAllTraces")
#   time_axis: Time samples (list)
#   n_traces: Number of traces 
#   """

#   # Plot kurtosis heatmap
#   plt.figure(figsize=figsize)
#   plt.imshow(kurt, aspect='auto', extent=(time_axis[0], time_axis[-1], n_traces-1, 0), 
#             cmap=cmap, vmin=vmin, vmax=vmax)
#   plt.colorbar()  

#   plt.title("Kurtosis")
#   plt.xlabel("Time [sec]")
#   plt.ylabel("Trace")  

#   # if plot_arrival==True:
#   #   # Plot arrivals
#   #   plt.scatter(np.array(tp), np.arange(n_traces), marker='|', s=1, color='red')
#   #   plt.scatter(np.array(ts), np.arange(n_traces), marker='|', s=1, color='red') 

def pickAllTraces(event, cut_trace, win, lpf=None, window=100,
                  save_file=None):
  """
  PS picking of all traces in the event data
  INPUT:
  event: Event data (TDMS object)
  cut_trace: Cut trace window
  win: Searching window for time arrivals. It is a LIST of TUPLES, for instance
    [(x0,y0), (x1,y1), ..., (xn, yn)] where n is the number of arrivals you want
    to find, xn is the the starting range, and yn is the ending range. In case
    you want to find P and S arrival, then your win will be [(xp,yp), (xs, ys)]
  lpf: Low-pass filter. 
    * Default is None: No filter is used
    * If used, specify as dictionary of f (low pass frequency) and fs (sampling
      frequency)
  window: Kurtosis calculation window. Default is 100.
  save_file: Saving arrivals to CSV file
    * Default is None: No file is saved
    * If saved, specify the filename to save
  
  OUTPUT:
  
  t: Autopicked arrival times 
  A: Autopicked arrival kurtosis

  NOTE: Both t and A are LIST, with size that depends on your specified "win"
        For two arrivals i.e. P and S arrivals, LIST will have shape (2,)
  """
  n_traces = len(event.zz)

  t, A = [], []
  for i in range(n_traces):
    x = kurtosisFindArrival(event=event, no_trace=i, cut_trace=cut_trace, 
                            win=win, lpf=lpf, window=window)
    t_, A_ = x
    t.append(t_)
    A.append(A_)

    if i%10==0:
      print("Finished picking until trace {}".format(i))

  if save_file!=None:
    # Save result to file
    df = pd.DataFrame({'D': event.zz})

    t, A = np.array(t), np.array(A)
    m, n = t.shape

    for i in range(n):
      # Record arrival times
      colname_t = 't'+'{}'.format(i+1)
      df[colname_t] = t[:,i]
      # Record arrival kurtosis
      colname_A = 'A'+'{}'.format(i+1)
      df[colname_A] = A[:,i]            

    df.to_csv(save_file, index=False)  
    print("Saved file to {}".format(save_file))
  
  return t, A
