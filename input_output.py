import numpy as np
import pandas as pd
import os
from datetime import datetime, timedelta
import glob
import pickle

def savePkl(data, filename):
  """
  Save data to Pickle

  INPUT:

  data: Numpy array
  filename: Name of file to be saved (extension .pkl)

  OUTPUT:

  Pickle file
  """
  with open(filename,'wb') as f:
    pickle.dump(data, f)

def loadPkl(filepath):
  """
  Read Pickle file
  
  INPUT:
  
  filepath: Path to file
  
  OUTPUT:
  
  Numpy array
  """
  with open(filepath, 'rb') as f:
    return pickle.load(f)

def readJMA_csv(filepath):
  """
  Read events from Japanese Meteorological Agency (JMA) catalog

  NOTE: Catalog must be in CSV file. Originally, it's a TXT. 
        Manually load it as CSV with Ms Excel. It has 11 columns: 
        'Date', 'Time', 'OTerr', 'Lat', 'LatErr', 'Long', 'LonErr', 'Dep', 
        'DepErr', 'Mag', 'Region'.
  
  INPUT:

  filepath: Path to JMA CSV catalog file

  OUTPUT:

  Pandas dataframe of catalogued events, columns are: 'index', 'Date', 'Time',
  'Lat', 'Lon', 'Dep', 'Mag', 'Loc'
  """
  df = pd.read_csv(filepath, header=0, dtype=str, index_col=False,
                  names=['Date', 'Time', 'OTerr', 'Lat', 'LatErr', 'Lon', 'LonErr',
                         'Dep', 'DepErr', 'Mag', 'Loc'])
  df = df.fillna('nan')

  # Convert to Pandas datetime
  df['Date'] = pd.to_datetime(df['Date'], format='%d/%m/%Y')
  df['Time'] = pd.to_datetime(df['Time'], format='%H:%M:%S.%f').dt.time

  # Edit magnitude, delete 'v'
  mag = [df.Mag.values[i][:3] for i in range(len(df))]
  mag = [np.nan if x=='na' else x for x in mag]
  
  mags = []
  for i in mag:
    if len(i)==3:
      # It's a magnitude, convert to float
      x = np.float(i)
    else:
      # It has unusual information e.g. 'na', convert to nan
      x = np.nan
    mags.append(x)

  df['Mag'] = mags

  # Convert other columns to float
  df['Lat'] = [float(df.Lat.values[i]) for i in range(len(df))]
  df['Lon'] = [float(df.Lon.values[i]) for i in range(len(df))]
  df['Dep'] = [float(df.Dep.values[i]) for i in range(len(df))]

  # Drop all NaN values of Lat and Lon
  catalog_df = df[['Date', 'Time', 'Lat', 'Lon', 'Dep', 'Mag', 'Loc']].\
               dropna(subset=['Lat']).reset_index(drop=True)
  return catalog_df

def getInfoFromJMA(filepath, df, utc=9, print_info=True):
  """
  Get information of a TDMS event file from a JMA catalog

  INPUT:

  filepath: Path to event file. The file must be in TDMS and have the following
    structure "connected whole_UTC_210501_153000.000" for event that happened in
    1 May 2021 at 15:30:00 UTC
  catalog_csv_path: Path to JMA catalog CSV file. JMA has the following key
    columns: "Date" and "Time"
  utc: UTC conversion. Default is 9 (UTC+9) for Japan
  print_info: Option to print info from the columns of catalog. Default is True.
    If False, it will return a dataframe. 
  """
  files = os.path.splitext(os.path.basename(filepath))[0]
  
  # Convert pandas datetime to string
  catalog_date, catalog_time = df.Date, df.Time
  catalog_date = catalog_date.apply(lambda x: x.strftime('%d/%m/%Y')).values
  catalog_time = catalog_time.apply(lambda x: x.strftime('%H:%M:%S')).values

  catalog_dt = [catalog_date[i]+' '+catalog_time[i][:5]+':00' for i in range(len(catalog_time))] 
  df['TDMSDatetime'] = catalog_dt # Add new column with catalog_dt

  # Extract timestamp from filename string
  timestamp = files[20:] # Omit connected whole bla bla ...
  
  # Convert string to datetime object
  timestamp = datetime.strptime(timestamp, '%Y%m%d_%H%M%S.%f')

  # Convert from UTC to local time
  timestamp = timestamp + timedelta(hours=utc)

  # Convert datetime object back to string
  timestamp = timestamp.strftime("%d/%m/%Y %H:%M:%S")

  # Find in catalog
  try:
    assert df.TDMSDatetime.str.contains(timestamp).any(), "no file"
    df = df[df.TDMSDatetime==timestamp]  
    if print_info==True:
      for i in range(len(df)):
        print('Info for file {}'.format(files))        
        print('Date           : {}'.format(df.Date.values[i]))
        print('Time           : {}'.format(df.Time.values[i]))
        print('Magnitude      : {}'.format(df.Mag.values[i]))
        print('Lat, Lon, Depth: ({}, {}, {})'.format(df.Lat.values[i], 
                                                     df.Lon.values[i],
                                                     df.Dep.values[i]))
        print('Location       : {}'.format(df.Loc.values[i]))
    if print_info==False:
      return df.drop(columns=['TDMSDatetime'])
  except:
    print('No info for file {}. Check in catalog.'.format(files))
    # return None
