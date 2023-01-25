# DosePainting

**Author**: Caterina Brighi

Scripts for the use of MRI data to derive heterogeneous dose prescriptions within the radiotherapy target, using patient-specific TCP data.

**Setup/Build/Install** 
To be able to use the above scipts you need to have the following python packages intalled:

```
- matplotlib
- SimpleITK
- numpy
- pandas
- datetime
- os
- glob
- gzip
- shutil
- xlsxwriter
- scipy
- dicom2nifti
- skimage
- pydicom
- platipy
```

You also need to have the *ImageAnalysisFunctions.py*, *ImageResultsVisualizationFunctions.py* and *ImageStatisticsFunctions.py* files in the same folders as the other .py files of this repo. 

**Usage** 
1. *dicom2nii.py* - this scripts convert the patients dicom files (MRI, RTSTRUCT, RTDOSE, CT) into nifti files.
2. *RegistrationCT_T2_DWI.py* - this script is used to register CT images to T2 MRI image space, and then to DWI image space, by means of a rigid registration using SimpleITK.
3. *alpha_beta_Protons.py* - this script is used to derive alpha and beta maps from maps of LET and PHYS dose, by applying the model developed in the publication [McNamara AL, Schuemann J, Paganetti H. A phenomenological relative biological effectiveness (RBE) model for proton therapy 
based on all published in vitro cell survival data. Phys Med Biol [Internet]. 2015 Nov 7;60(21):8399–416.](https://iopscience.iop.org/article/10.1088/0031-9155/60/21/8399)).
4. *Dpresc_fromTCP.py* - this script uses patient-specific information of tumour cellularity (derived from DWI-ADC data), alpha and beta to derive an optimal dose prescription based on a tumour control probability model (based on poissonian distribution and LQ model).
5. *CTV_doseStats.py* - this script generates an excel spreadsheet with some statistics metrics calculated on the PHYS dose maps and the optimised dose prescriptions in the clinical target volume (CTV). These metrics include: Mean, std, median, max and min dose.

**Directory Structure** 
The original data directory before converting the dicom images to nifti with the *dicom2nii.py* script should be the following:
```
Data
 dicom
   MRI
     Patient ID
       Study date
         MR
           Sequence X
             file1.dcm
             file2.dcm
             ...
       ...
     ...
   RTDATA
     Patient ID
       CT
         file1.dcm
         file2.dcm
         ...
       LEM
         file1.dcm
         file2.dcm
         ...
       LET
         file1.dcm
         file2.dcm
         ...
       MKM
         file1.dcm
         file2.dcm
         ...
       PHYS
         file1.dcm
         file2.dcm
         ...
       PLAN
         file1.dcm
         file2.dcm
         ...
       RTSTRUCT
         file1.dcm
         file2.dcm
         ...
 ```       

**Citation**
If you use code from this repository, please cite [to be added]
