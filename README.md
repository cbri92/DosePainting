# DosePainting

**Author**: Caterina Brighi

Scripts for the use of MRI data to derive heterogeneous dose prescriptions within the radiotherapy target, using patient-specific TCP data.

**Setup/Build/Install** 
To be able to use the above scipts you need to have the following python packages intalled:

-matplotlib
-SimpleITK
-numpy
-pandas
-datetime
-os
-glob
-gzip
-shutil
-xlsxwriter
-scipy
-dicom2nifti
-skimage
-pydicom
-platipy

You also need to have the *ImageAnalysisFunctions.py*, *ImageResultsVisualizationFunctions.py* and *ImageStatisticsFunctions.py* files in the same folders as the other .py files of this repo. 

**Usage** 
-*dicom2nii.py* - this scripts convert the patients dicom files (MRI, RTSTRUCT, RTDOSE, CT) into nifti files.
-*alpha_beta_Protons.py* - this script is used to derive alpha and beta maps from maps of LET and PHYS dose, by applying the model developed by McNamara et al. 


**Directory Structure** 
NA

**Citation**
If you use code from this repository, please cite [to be added]
