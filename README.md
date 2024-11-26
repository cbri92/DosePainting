# DosePainting

**Author**: Caterina Brighi

Scripts for the use of MRI and RT data to derive heterogeneous dose prescriptions within the radiotherapy target, using patient-specific cellularity data.

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

IMPORTANT: You always need to have the *ImageAnalysisFunctions.py* and *ImageStatisticsFunctions.py* files, which you can find in the utils folder, in the same folders as the other .py files when running pre-processing and analysis steps (i.e. Pre-processing and Analysis folders). 

**Directory Structure** 
The original data directory before running the pre-processing steps should be the following:
```
Data supradirectory
  Patient directory
    dicom
      CT
        file1.dcm
        file2.dcm
        ...
      RTDOSE
        RBE_TOT.dcm
      RTSTRUCT
        rtss_xxxx.dcm
       
    MRI
      baseline
        micro
          ADC_50_400_1000.nii
          cellApp_CORRECTED_CTV.nii
          cellApp_CORRECTED_GTV.nii
          gtv_voxel_ok_3D.nii
          ctv_voxel_ok_3D.nii
  
   regFiles
      rigid_CTonT2.txt
      rigid_T2onDWI.txt
      xf_CTonT2.txt
      xf_T2onDWI.txt
 ```       

NOTE: The nifti files in the MRI/baseline/micro folder are calculated elsewhere by applying a model to convert DWI MRI ADC data into cellularity data (see publication: ). The text files in the regFiles folder are also generated in ITK-SNAP by performing a manual/semiautomated rigid registration of the CT to T2 image and of the T2 image to the DWI image.

**Citation**
If you use code from this repository, please cite [to be added]
