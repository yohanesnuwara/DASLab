# DASLab

**What the heck is DAS?** DAS, stands for Distributed Acoustic Sensing, is a geophysical measurement technology that uses laser light scattering inside fiber optic cables to measure dynamic strain and seismic waves. The strain measured by DAS is analogous to pressure measured by seismic geophones. Unlike geophones which have spacing between each, DAS uses the whole fiber. This is the reason why DAS is used in seismic acquisition. Many believe, like [Bona and Pevzner (2018)](https://www.tandfonline.com/doi/abs/10.1071/ASEG2018abW8_4F), that DAS has many advantages such as consistent amplitudes, less tube wave, clearer reflections, faster and cheaper cable acquisition.  

## Check out the notebooks to see how we make use of DASLab

* **DASAutopick**: Displaying DAS data and using multiple algorithms such as Kurtosis, AIC/BIC, STA/LTA, and Short Time Fourier Transform (STFT) for automated detection and picking of P-S wave events from DAS recordings. 
* **DASCompare**: Using waveform visualization, STFT spectogram, and f-x waterfall plot to compare the difference in straight and helical DAS cable recordings.

## Useful bags

* [Global subsea fiber optic network map](https://submarine-cable-map-2019.telegeography.com/)
* [Geothermal data @ GDR OpenEi](https://gdr.openei.org/submissions/980)
* [DAS technical explanations](http://docs.energistics.org/PRODML/PRODML_TOPICS/PRO-DAS-000-016-0-C-sv2000.html)
* [Seismo-Live](http://seismo-live.org/)
* [Eileenmartin](https://github.com/eileenrmartin) lots of repo about DAS
* Distpy (Schlumberger)
* Ariel Lellouch, [DAS data](https://github.com/ariellellouch/DASDetection)
* [DAS anomaly detection ML](https://github.com/rroy1212/DAS_Anomaly_Detection)

## Dependencies:
* [npTDMS](https://pypi.org/project/npTDMS/0.25.0/): For reading TDMS files (**required version 0.25.0**)
* [ObsPy](https://pypi.org/project/obspy/): For seismological operations 
* [Utm](https://pypi.org/project/utm/): For converting lat-long to UTM coordinates
