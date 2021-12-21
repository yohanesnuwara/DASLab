# DASLab

**What the heck is DAS?** DAS, stands for Distributed Acoustic Sensing, is an emerging geophysical technology that uses laser light scattering inside fiber-optic cables to measure dynamic strain and seismic waves. The strain measured by DAS is analogous to pressure measured by seismic geophones. Unlike geophones which have spacing between each, DAS uses the whole fiber. This is the reason why DAS is applied in many monitoring and profiling applications, such as earthquake (seismology), seismic data acquisition and reservoir monitoring of oil fields, construction, wind turbine, chemical tank storage, dam, and highways. Many believe, like [Bona and Pevzner (2018)](https://www.tandfonline.com/doi/abs/10.1071/ASEG2018abW8_4F), that DAS has many advantages such as consistent amplitudes, less tube wave, clearer reflections, faster and cheaper cable acquisition.  

**DASLab** is a repository that contains programs, modules, and functions for analysis and processing of DAS fiber-optic data that I wrote during my research with the [CO2 Storage Research Group at RITE](http://www.rite.or.jp/co2storage/en/) in Japan. 

## Check out the notebooks to see how we make use of DASLab

* **DASAutopick**: Displaying DAS data, Butterworth bandpass filtering, and using multiple algorithms such as Kurtosis, AIC/BIC, STA/LTA, and Short Time Fourier Transform (STFT) for automated detection and picking of P-S wave events from DAS recordings. 
* **DASCompare**: Using waveform visualization, STFT spectogram, and f-x waterfall plot to compare the difference in straight and helical DAS cable recordings.
* **DASContinuous**: Applying the automated picking algorithms on a simulated continuous DAS recordings.
* **DASDiagnostics**: Applying trace normalization and spectral analysis using f-x waterfall plot and F-K filter plot.
* **DASCatalogMap**: Retrieving information from JMA catalog into Pandas DataFrame and plotting events on a map and polar plot.
* **DASDetectivity**: Analyzing the detectivity of DAS to natural earthquakes based on picked arrival times, amplitude or Kurtosis plot versus epicentral distance, and azimuthal dependence on detectivity on polar plot. 
* **DASSignalNoise**: Calculation of background noise level, signal-to-noise ratio (SNR) analysis, and to explore trace operations such as stacking and integration to improve SNR.

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
