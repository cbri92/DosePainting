# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 09:50:05 2023

@author: Caterina Brighi

This script:
    1. applies the inverse transfor to bring the dose painting prescriptions from DWI to CT space
    2. assigns values of 74 Gy dose to CTV margins (CTV-GTV)
    3. converts the dose painting prescriptions into inverse dose painting prescriptions
    4. converts the inverse dose paitning prescriptions from nifti to dicom RT dose files
"""


#%% Import functions 

import SimpleITK as sitk
import os
import glob
from ImageAnalysisFunctions import getNonZeroStats, set_mask_value, generate_mask, flip_image
from ConvertNii_ToDoseFiles import * #Use this functions only when using the dicom CT series file as reference dicom

#%% Set Working directory
        
data_supradir = 'path/to/patients/data/supra/directory/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    print('Creating inverse dose prescriptions for '+current)
    subj_dir = data_supradir+current
    subj_name = current
    
    #Read directories' path
    dose_dir = data_supradir+current+'/MRI/baseline/Prescribed_dose'
    reg_dir = data_supradir+current+'/regFiles'
    
    #Make a directory for new inverse dose images
    if not os.path.exists(data_supradir+current+'/RTPRESC'):#if it does not already exist, create a directory where the dose prescription files will be saved
        os.mkdir(data_supradir+current+'/RTPRESC')
    prsc_dir = data_supradir+current+'/RTPRESC'
    
    #Read CT image
    CT = sitk.ReadImage(subj_dir+'/ct.nii')
    
    #Read dose prescription image
    DP_inDWI = sitk.ReadImage(dose_dir +'/Dose_painted_GTV_final.nii')
    
    #%%Read transformation files
    transform1 = sitk.ReadTransform(reg_dir+'/xf_CTonT2.txt')
    versor1 = sitk.VersorRigid3DTransform(transform1)
    
    transform2 = sitk.ReadTransform(reg_dir+'/xf_T2onDWI.txt')
    versor2 = sitk.VersorRigid3DTransform(transform2)

    composite_transform = sitk.CompositeTransform([versor2, versor1])
    inverse_transform = composite_transform.GetInverse()
    
    #%%Apply inverse transform to Dose prescriptions to bring them in CT space from DWI space
    DP_inCT = sitk.Resample(DP_inDWI, CT, inverse_transform, sitk.sitkNearestNeighbor)
   
    #%%Calculate Maximum dose prescription values    
    Stats_dose = getNonZeroStats(DP_inCT)
    MAX_dose = int(Stats_dose['Max intensity'])
    
    #%%Generate inverse dose prescription as Dinv = Dmax - Dpresc, where Dmax = maximum dose estimated for each patient
    Inv_DP_inCT = MAX_dose-DP_inCT
    
    #%%Read CTV HD
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'CTV*'):
        if (('CTV_7' in filename) or ('CTV7' in filename) or ('CTV_high' in filename) or ('CTV_HD' in filename)):
            CTV_HD_path = filename
                        
    CTV_HD = sitk.ReadImage(CTV_HD_path) #read PTV HD
  
    #%%Generate HD margins
    Valid_presc = DP_inCT>0
    HD_margin = (CTV_HD-Valid_presc)>0
    sitk.WriteImage(HD_margin, prsc_dir+'/HD_margin.nii')
    
    #%%Assign 74 Gy to all voxels in HD margins
    DP_inCT_inCTV = set_mask_value(DP_inCT, HD_margin, 74)
    
    Inv_DP_inCT_inCTV = set_mask_value(Inv_DP_inCT, HD_margin, (MAX_dose-74))
    Inv_DP_inCT_inCTV = generate_mask(Inv_DP_inCT_inCTV, CTV_HD)
    
    #%%Write dose prescriptions and inverse dose prescriptions        
    sitk.WriteImage(DP_inCT_inCTV, prsc_dir +'/Dpainted.nii')
    sitk.WriteImage(Inv_DP_inCT_inCTV, prsc_dir +'/Inv_Dpainted.nii')
    
    #%%Convert inverse dose prescriptions into dcm RT dose files using CT dicom series as reference dicom
    
    #First we need to flip the inverse prescriptins wrt y axis, as the dii2nii algo used originally to convert all the images flipped this wrt to original dicom images
    Inv_DP_inCT_inCTV = flip_image(Inv_DP_inCT_inCTV, 'y')
    ct_ref = subj_dir+'/dicom/CT/' #Path to directory of CT dicom series
    
    convert_nii_to_dicom_RTdosefile(Inv_DP_inCT_inCTV, ct_ref, output_directory=prsc_dir, out_filename="Inv_Dpainted.dcm")
