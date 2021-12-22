<p align="left">
  <img src="https://user-images.githubusercontent.com/51282928/146875829-05c3c234-35df-4646-8d17-1435dab6996e.png" width="200" />
</p>

[![DOI](https://zenodo.org/badge/365670829.svg)](https://zenodo.org/badge/latestdoi/365670829)

**What is DAS?** DAS, stands for Distributed Acoustic Sensing, is an emerging geophysical technology that uses laser light scattering inside fiber-optic cables to measure dynamic strain and seismic waves. To learn more detailed about how the sensing works, I recommend this [blog](https://motionsignaltechnologies.com/what-is-das-and-what-is-it-measuring/). DAS is widely applied in many monitoring and profiling applications, such as earthquake (seismology), seismic data acquisition and reservoir monitoring of oil fields, construction, wind turbine, chemical tank storage, dam, and highways. Here, I wrote a comprehensive [article](https://www.linkedin.com/pulse/distributed-acoustic-sensing-answer-sustainability-quest-nuwara/) about DAS. Many believe, like [Bona and Pevzner (2018)](https://www.tandfonline.com/doi/abs/10.1071/ASEG2018abW8_4F), that DAS has many advantages such as consistent amplitudes, less tube wave, clearer reflections, faster and cheaper cable acquisition.  

<p align="center">
  <img src="https://user-images.githubusercontent.com/51282928/146873724-264c0aa5-ba7f-41ae-87c2-f88a03c534ba.png" width="400" />
</p>

**DASLab** is a repository that contains programs, modules, and functions for analysis and processing of DAS fiber-optic data that I developed during my research with the [CO2 Storage Research Group at RITE](http://www.rite.or.jp/co2storage/en/) in Japan. 

## Check out the notebooks to see how we make use of DASLab

* **[DASAutopick](https://github.com/yohanesnuwara/DASLab/blob/main/notebooks/DASAutopick.ipynb)**: Displaying DAS data, Butterworth bandpass filtering, and using multiple algorithms such as Kurtosis, AIC/BIC, STA/LTA, and Short Time Fourier Transform (STFT) for automated detection and picking of P-S wave events from DAS recordings. 
* **[DASCompare](https://github.com/yohanesnuwara/DASLab/blob/main/notebooks/DASCompare.ipynb)**: Using waveform visualization, STFT spectogram, and f-x waterfall plot to compare the difference in straight and helical DAS cable recordings.
* **[DASContinuous](https://github.com/yohanesnuwara/DASLab/blob/main/notebooks/DASContinuous.ipynb)**: Applying the automated picking algorithms on a simulated continuous DAS recordings.
* **[DASSignalNoise](https://github.com/yohanesnuwara/DASLab/blob/main/notebooks/DASSignalNoise.ipynb)**: Calculation of background noise level, signal-to-noise ratio (SNR) analysis, and to explore trace operations such as stacking and integration to improve SNR.
* **[DASDiagnostics](https://github.com/yohanesnuwara/DASLab/blob/main/notebooks/DASDiagnostics.ipynb)**: Applying trace normalization and spectral analysis using f-x waterfall plot and F-K filter plot.
* **[DASCatalogMap](https://github.com/yohanesnuwara/DASLab/blob/main/notebooks/DASCatalogMap.ipynb)**: Retrieving information from JMA catalog into Pandas DataFrame and plotting events on a map and polar plot.
* **[DASDetectivity](https://github.com/yohanesnuwara/DASLab/blob/main/notebooks/DASDetectivity.ipynb)**: Analyzing the detectivity of DAS to natural earthquakes based on picked arrival times, amplitude or Kurtosis plot versus epicentral distance, and azimuthal dependence on detectivity on polar plot. 

## Cite this repository

If you use this repository, please cite it as below:

> Nuwara, Y. (2021). DASLab: Distributed Acoustic Sensing Lab (Version 1.0.0) [Computer software]

## Dependencies:
* [npTDMS](https://pypi.org/project/npTDMS/0.25.0/): For reading TDMS files (**required version 0.25.0**)
* [ObsPy](https://pypi.org/project/obspy/): For seismological operations 
* [Utm](https://pypi.org/project/utm/): For converting lat-long to UTM coordinates

<!--
## Useful bags

* [Global subsea fiber optic network map](https://submarine-cable-map-2019.telegeography.com/)
* [Geothermal data @ GDR OpenEi](https://gdr.openei.org/submissions/980)
* [DAS technical explanations](http://docs.energistics.org/PRODML/PRODML_TOPICS/PRO-DAS-000-016-0-C-sv2000.html)
* [Seismo-Live](http://seismo-live.org/)
* [Eileenmartin](https://github.com/eileenrmartin) lots of repo about DAS
* Distpy (Schlumberger)
* Ariel Lellouch, [DAS data](https://github.com/ariellellouch/DASDetection)
* [DAS anomaly detection ML](https://github.com/rroy1212/DAS_Anomaly_Detection)
