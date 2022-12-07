# MISARA:Matlab Interface for Seismo-Acoustic aRray Analysis
MISARA (Matlab Interface for Seismo-Acoustic aRray Analysis) is a open-source Matlab-based GUI designed to perform analyses of seismic and acoustic waveform data. A suite of well-established algorithms for volcano seismic and acoustic signal processing have been integrated into our GUI interface, with a special focus on array techniques (for more details, see Rost and Thomas, 2000). We note that although MISARA was developed to facilitate the analysis of seismic and acoustic signals in volcanic environments, it can be used for other research purposes. Furthermore, owing to its modular structure, it is possible to easily integrate additional functionalities.
The different data analysis modules of MISARA are independent of each other. The modules were designed to easily manage every step of the data processing and to quickly inspect the results. Most of the processes are automated, reducing user’s errors and efforts. One advantage consists of the possibility to reset some parameters directly from the module itself, allowing to repeat the analysis many times. Other fundamental aspects of this modular structure are the possibility to deal with different formats of input traces, the systematic saving of the results and the optional activation of many subroutines.
The main structure of the interface consists of:

-	Home window, the main panel for the management of all utilities of MISARA.
-	Data preparation window, for the formatting of the Input data.
-	Data Pre-processing modules, for the data quality control.
-	Signal Features modules, for those analytic routines that support the array methods, such as spectral, amplitude, polarization and detection analysis.
-	Array analysis modules, for the source localization methods based on the multichannel techniques.

# Requirements
MISARA can be run on any operation system with Matlab from Release 2021b. In addition, you need to have the following toolboxes installed: 

-	Control System Toolbox, version 10.5.
-	Financial Toolbox, version 5.12.
-	Mapping Toolbox, version 4.7.
-	Signal Processing Toolbox, version 8.1.
-	Statistics and Machine Learning Toolbox, version 11.4.
-	Wavelet Toolbox, version 5.1.

# Documentations
You can consult the manual of the software from the "Help" menu of the control panel of MISARA or you can open it from "MISARA/Doc/" directory.

# Starting MISARA
After you have downloaded the software, you should unzip the source code to a suitable directory. To start MISARA, you should run “MISARA.m” by pressing F5 in the Matlab editor. The software automatically sets the main paths and functions, but it requires the installation of the “GIPPtools” and “irisFetch” libraries. For any information, see the section 2.1 of the User Manual. After the running, the Home Window appears on the screen, allowing you to set the analysis parameters and to use all modules. To test the software, you can refer to “Video Tutorial” from the “Help” menu located in the upper-right part of the Home Window. To use these video tutorials, we suggest you to download the software from this URL: https://doi.org/10.5281/zenodo.7410076. 

The purpose of these video tutorials is to train users to perform different types of analyses of seismic and acoustic waveform data acquired in volcanic environment. In particular they are grouped into three sections. They show you a series of brief tutorials on how to use the software on three real cases studies. In First one, we will perform the analysis of volcanic tremor recorded by a seismic array deployed at Mt. Etna (Italy) in 2011, when the volcano produced intense lava fountain activity from its New South East Crater (NSEC; for more details about volcanic activity, see Bencke et al, 2014). In the second one, we will demonstrate analyses of Long Period (LP) and Very Long Period (VLP) earthquakes recorded by Mt. Etna permanent seismic network in 2010, accompanying explosive activity at the Bocca Nuova crater (BN; for more details about volcanic activity, see Andronico et al, 2010). In the third one, we will show how to analyse the infrasound data acquired by an infrasound array deployed at Mt. Etna in 2019, when the NSEC crater was affected by intense Strombolian activity (for more details about volcanic activity, see De Angelis et al, 2020). The raw data that we will use in these tutorials are in the “MISARA/Data_org/” directory. However, it is possible to use the Matlab data format located in “MISARA/Data_example/” directory. To use this demostration dataset, we suggest you to download the software from this URL: https://doi.org/10.5281/zenodo.7410076.  

# Citation 
This is a modified version of the GSpecDisp package (Sadeghisorkhani et al., 2017).
If you use this code for your work, please cite the following DOI:
-	https://doi.org/10.5281/zenodo.7410076

# Contact
You can send an email to vittorio.minio@phd.unict.it to report suggestions, comments and bugs.

# References
-	Andronico, D., Lo Castro, M. D., Sciotto, M., Spina, L., 2013. The 2010 ash emissions at the summit craters of Mt Etna: Relationship with seismo‐acoustic signals. Journal of Geophysical Research: Solid Earth, 118(1), 51-70.
-	Behncke, B., Branca, S., Corsaro, R. A., De Beni, E., Miraglia, L., Proietti, C., 2014. The 2011–2012 summit activity of Mount Etna: Birth, growth and products of the new SE crater. Journal of Volcanology and Geothermal Research, 270, 10-21.
-	De Angelis, S., Haney, M. M., Lyons, J. J., Wech, A., Fee, D., Diaz-Moreno, A., Zuccarello, L., 2020. Uncertainty in detection of volcanic activity using infrasound arrays: examples from Mt. Etna, Italy, Frontiers in Earth Science, 8, 169.
-	Rost, S., and Thomas, C., 2002. Array seismology: Methods and applications. Reviews of geophysics, 40(3).
-	Sadeghisorkhani, H., Gudmundsson, O., Tryggvason, A., 2017. GSpecDisp: a Matlab GUI package for phase-velocity dispersion measurements from ambient-noise correlations, Computers and Geosciences.
