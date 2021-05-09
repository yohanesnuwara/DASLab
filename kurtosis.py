import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

def kurtosis_find_arrival(event, no_trace, cut_trace, win1, win2, lpf=None,
                          window=100, plot=False, title=None):
  """
  Find P-S arrival using Kurtosis and plot

  INPUT:

  event: Event data (TDMS object)
  cut_trace: Cut trace window. Specify as tuple (start,end) with start is the
    start time and end is the end time
  win1: Searching window for P arrivals. Specify as tuple (start,end) with start is the
    start time and end is the end time
  win2: Searching window for S arrivals. Specify as tuple (start,end) with start is the
    start time and end is the end time
  lpf: Low-pass filter. 
    * Default is None: No filter is used
    * If used, specify as dictionary of f (low pass frequency) and fs (sampling
      frequency)
  window: Kurtosis calculation window. Default is 100.
  plot: Option to plot the result. Default is False, so no plot is produced, but
    the following OUTPUTS will be produced.

  OUTPUT:

  tp, ts: P and S arrival of all traces (1D arrays)
  Ap, As: Kurtosis value when P and S arrivals of all traces (1D arrays)
  """

  # Cut trace 
  tcut, ycut = cutTrace(event, no_trace, cut_trace)
  # t, z = event.tt, event.zz # axes
  # index_cut = [int(np.where(t==cut_trace[0])[0]), int(np.where(t==cut_trace[1])[0])]
  # tcut = t[index_cut[0]:index_cut[1]]

  # fibre_data = event.data
  # y = fibre_data[:,no_trace] 
  # ycut = y[index_cut[0]:index_cut[1]]

  if lpf!=None:
    # Low pass filter is applied
    f, fs = lpf["f"], lpf["fs"]
    ycut = butter_lowpass_filter(ycut, f, fs, order=5)

  """
  KURTOSIS
  """

  # def kurtosis(data, window=10):
  #   """ Calculate kurtosis over a specified window """
  #   kurt = []
  #   for i in range(len(data)):
  #     a = data[i:i+window]
  #     std, mean = np.std(data), np.mean(data)
  #     y = np.sum((a - mean)**4) / window
  #     k = y / std**4
  #     kurt.append(k)
  #   return kurt

  # Kurtosis calculation
  ycut_kurt = kurtosis(ycut, window=window)

  """
  IDENTIFY P AND S
  """

  # Divide into 2 windows
  win1 = np.where((tcut >= win1[0]) & (tcut <= win1[1]))[0]
  win2 = np.where((tcut >= win2[0]) & (tcut <= win2[1]))[0]
  win1_start, win1_end = win1[0], win1[-1]
  win2_start, win2_end = win2[0], win2[-1]

  win1_trace, win2_trace = ycut_kurt[win1_start:win1_end], ycut_kurt[win2_start:win2_end]

  # Search the highest point
  max_win1, max_win2 = max(win1_trace), max(win2_trace)
  Ap, As = max_win1, max_win2

  # Search time at the highest point (arrivals)
  tp, ts = np.where(ycut_kurt == max_win1)[0], np.where(ycut_kurt == max_win2)[0]
  tp, ts = tcut[tp][0], tcut[ts][0]

  if plot==True:
    # Plot seismogram and kurtosis
    plt.figure(figsize=(11,8))

    plt.subplot(2,1,1)
    plt.plot(tcut, ycut, color="black")
    plt.title("Trace Cut No. {} of {}".format(no_trace, title))
    # plt.xlim(tcut[0], tcut[-1])
    plt.xlabel("Time [sec]")
    plt.ylabel("Amplitude")
    plt.grid()

    plt.axvline(tp, color="blue"); plt.axvline(ts, color="red")

    plt.subplot(2,1,2)
    plt.title("Kurtosis of Trace No. {}".format(no_trace))
    plt.plot(tcut, ycut_kurt, color="red")
    # plt.xlim(10,20)
    plt.xlabel("Time [sec]")
    plt.ylabel("Kurtosis")
    plt.grid()

    plt.tight_layout(1.1)
    plt.show()
  
  if plot==False:
    # No plot is produced
    return tp, ts, Ap, As

def kurtosisHeatmap(kurt, time_axis, n_traces, cmap='jet', figsize=(7,5),
                    vmin=None, vmax=None):
  """
  Plot heatmap of kurtosis of all traces

  INPUT:

  kurt: 2D array of calculated kurtosis (output of running "pickAllTraces")
  time_axis: Time samples (list)
  n_traces: Number of traces 
  """

  # Plot kurtosis heatmap
  plt.figure(figsize=figsize)
  plt.imshow(kurt, aspect='auto', extent=(tcut[0], tcut[-1], n_traces-1, 0), 
            cmap=cmap, vmin=vmin, vmax=vmax)
  plt.colorbar()  

  plt.title("Kurtosis")
  plt.xlabel("Time [sec]")
  plt.ylabel("Trace")  

  # if plot_arrival==True:
  #   # Plot arrivals
  #   plt.scatter(np.array(tp), np.arange(n_traces), marker='|', s=1, color='red')
  #   plt.scatter(np.array(ts), np.arange(n_traces), marker='|', s=1, color='red')  

def pickAllTraces(event, cut_trace, win1, win2, lpf=None, window=100,
                  save_file=None):
  """
  PS picking of all traces in the event data

  INPUT:

  event: Event data (TDMS object)
  cut_trace: Cut trace window
  win1: Searching window for P arrivals
  win2: Searching window for S arrivals
  lpf: Low-pass filter. 
    * Default is None: No filter is used
    * If used, specify as dictionary of f (low pass frequency) and fs (sampling
      frequency)
  window: Kurtosis calculation window. Default is 100.
  save_file: Saving arrivals to CSV file
    * Default is None: No file is saved
    * If saved, specify the filename to save
  
  OUTPUT:

  tp, ts: P and S arrival of all traces (1D arrays)
  Ap, As: Kurtosis value when P and S arrivals of all traces (1D arrays)
  """
  n_traces = len(event.zz)

  tp, ts, Ap, As = [], [], [], []
  for i in range(n_traces):
    x = kurtosis_find_arrival(event=event, no_trace=i, cut_trace=cut_trace, 
                              win1=win1, win2=win2, lpf=lpf, window=window)
    tp.append(x[0])
    ts.append(x[1])
    Ap.append(x[2])
    As.append(x[3])

    if i%10==0:
      print("Finished picking until trace {}".format(i))

  if save_file!=None:
    # Save result to file
    df = pd.DataFrame({'D': event.zz, 'tp': tp, 'ts': ts, 'Ap': Ap, 'As': As})
    df.to_csv(save_file, index=False)  
    print("Saved file to {}".format(save_file))
  
  return tp, ts, Ap, As
