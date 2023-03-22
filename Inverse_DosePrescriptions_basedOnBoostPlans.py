# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 09:50:05 2023

@author: Caterina Brighi

This script:
    1. applies the inverse transfor to bring the dose prescriptions from DWI to CT space
    2. assigns values of min dose to PTV margins (PTV- CTV)
    3. converts the dose prescriptions into inverse dose prescriptions 
"""


#%% Import functions 

import SimpleITK as sitk
import os
import glob
from ImageAnalysisFunctions import *
from ConvertNii_ToDoseFiles import * #Use this functions only when using the dicom CT series file as reference dicom

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

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
    DP_noBoneDWI = sitk.ReadImage(dose_dir +'/Dose_optimised_CTVnoBone.nii')
    
    #%%Read transformation files
    transform1 = sitk.ReadTransform(reg_dir+'/xf_CTonT2.txt')
    versor1 = sitk.VersorRigid3DTransform(transform1)
    
    transform2 = sitk.ReadTransform(reg_dir+'/xf_T2onDWI.txt')
    versor2 = sitk.VersorRigid3DTransform(transform2)

    composite_transform = sitk.CompositeTransform([versor2, versor1])
    inverse_transform = composite_transform.GetInverse()
    
    #%%Apply inverse transform to Dose prescriptions to bring them in CT space from DWI space
    DP_noBoneCT = sitk.Resample(DP_noBoneDWI, CT, inverse_transform, sitk.sitkNearestNeighbor)
   
    #%%Calculate Max and Min dose prescription values    
    Stats_dose_noBone = getNonZeroStats(DP_noBoneCT)
    MAX_dose_noBone = int(Stats_dose_noBone['Max intensity'])
    
    #%%Generate inverse dose prescription as Dinv = Dmax - Dpresc, where Dmax = maximum dose estimated for each patient
    Inv_DP_noBoneCT = MAX_dose_noBone-DP_noBoneCT
    
    #%%Read TV HD and PTV LD
    #Read PTV LD
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'PTV*'):
        if (('PTV_5' in filename) or ('PTV5' in filename) or ('PTV_3' in filename) or ('PTV_low' in filename) or ('PTV_LD' in filename)):
            PTV_LD_path = filename
                        
    PTV_LD = sitk.ReadImage(PTV_LD_path) #read PTV LD
    
    #Read PTV HD
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'PTV*'):
        if (('PTV_7' in filename) or ('PTV7' in filename) or ('PTV_high' in filename) or ('PTV_HD' in filename)):
            PTV_HD_path = filename
                        
    PTV_HD = sitk.ReadImage(PTV_HD_path) #read PTV HD
    
    #%%Generate HD and LD margins
    
    Valid_presc = DP_noBoneCT>0
    HD_margin = (PTV_HD-Valid_presc)>0
    LD_margin = (PTV_LD-PTV_HD)
    sitk.WriteImage(HD_margin, prsc_dir+'/HD_margin.nii')
    sitk.WriteImage(LD_margin, prsc_dir+'/LD_margin.nii')
    
    #%%Generate mask of dose prescriptions and inverse dose prescriptions in CTV    
    DP_noBoneCT_inCTV = generate_mask(DP_noBoneCT, Valid_presc)    
    Inv_DP_noBoneCT_inCTV = generate_mask(Inv_DP_noBoneCT, Valid_presc)
    
    #%%Assign 74 Gy to all voxels in HD_margins, and 54 Gy to all voxels in LD_margins
    DP_noBoneCT_inPTV = set_mask_value(DP_noBoneCT_inCTV, HD_margin, 74)  
    DP_noBoneCT_inPTV = set_mask_value(DP_noBoneCT_inPTV, LD_margin, 54)
    
    Inv_DP_noBoneCT_inPTV = set_mask_value(Inv_DP_noBoneCT_inCTV, HD_margin, (MAX_dose_noBone-74))
    Inv_DP_noBoneCT_inPTV = set_mask_value(Inv_DP_noBoneCT_inPTV, LD_margin, (MAX_dose_noBone-54))
    
    #%%Write dose prescriptions and inverse dose prescriptions        
    sitk.WriteImage(DP_noBoneCT_inPTV, prsc_dir +'/DP_noBone_CfBoost.nii')
    sitk.WriteImage(Inv_DP_noBoneCT_inPTV, prsc_dir +'/Inv_DP_noBone_CfBoost.nii')
    
    #%%Convert inverse dose prescriptions into dcm RT dose files using CT dicom series as reference dicom
    
    #First we need to flip the inverse prescriptins wrt y axis, as the dii2nii algo used originally to convert all the images flipped this wrt to original dicom images
    Inv_DP_noBoneCT_inPTV = flip_image(Inv_DP_noBoneCT_inPTV, 'y')
    ct_ref = subj_dir+'/CT/' #Path to directory of CT dicom series
    
    convert_nii_to_dicom_RTdosefile(Inv_DP_noBoneCT_inPTV, ct_ref, output_directory=prsc_dir, out_filename="Inv_DP_noBone_CfBoost.dcm")
