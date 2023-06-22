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
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti/' #Set working directory
dcm_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/dicom/RTDATA/'

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    print('Creating inverse dose prescriptions for '+current)
    subj_dir = data_supradir+current
    subj_name = current
    
    MRI_dir = subj_dir+'/MRI/'

    studies = [ s.name for s in os.scandir(MRI_dir) if s.is_dir() ]
    studies = [x.replace('study_', '') for x in studies] #Remove word "study_" from list of studies
    studies2 = [int(x) for x in studies] #Convert list of strings in list of integers
    small_study = min(studies2) #Find the smallest study date
    study = 'study_'+str(small_study) #Define which one is the baseline (earliest) study
    
    #Read directories' path
    dose_dir = MRI_dir+study+'/Prescribed_dose'
    reg_dir = MRI_dir+study+'/regFiles'
    
    #Make a directory for new inverse dose images
    if not os.path.exists(data_supradir+current+'/RTPRESC'):#if it does not already exist, create a directory where the dose prescription files will be saved
        os.mkdir(data_supradir+current+'/RTPRESC')
    prsc_dir = data_supradir+current+'/RTPRESC'
    
    #Read CT image
    CT = sitk.ReadImage(subj_dir+'/ct.nii')
    
    #Read dose prescription image
    DP_noBoneDWI = sitk.ReadImage(dose_dir +'/Dose_optimised_CTVnoBone.nii')
    
    #%%Read transformation files
    transform1 = sitk.ReadTransform(reg_dir+'/xf_CTonT2_1.txt')  
    transform2 = sitk.ReadTransform(reg_dir+'/xf_CTonT2_2.txt')

    composite_transform = sitk.CompositeTransform([transform2, transform1])
    inverse_transform = composite_transform.GetInverse()
    
    #%%Apply inverse transform to Dose prescriptions to bring them in CT space from DWI space
    DP_noBoneCT = sitk.Resample(DP_noBoneDWI, CT, inverse_transform, sitk.sitkNearestNeighbor)
   
    #%%Calculate Max and Min dose prescription values    
    Stats_dose_noBone = getNonZeroStats(DP_noBoneCT)
    MAX_dose_noBone = int(round(Stats_dose_noBone['Max intensity'],1))
     
    #%%Generate inverse dose prescription as Dinv = Dmax - Dpresc, where Dmax = maximum dose estimated for each patient
    Inv_DP_noBoneCT = MAX_dose_noBone-DP_noBoneCT
    
    #%%Read TV HD and PTV LD
    #Read PTV LD
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'PTV*'):
        if (('PTV_4' in filename) or ('PTV_bassa' in filename) or ('PTV_39.6.' in filename) or ('PTV_39.6_4mm' in filename) or ('PTV3' in filename) or ('PTV_low' in filename) or ('PTV_LD.' in filename) or ('PTVLD.' in filename)):
            PTV_LD_path = filename
                        
    PTV_LD = sitk.ReadImage(PTV_LD_path) #read PTV LD
    
    #Read PTV HD
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'PTV*'):
        if (('PTV_70.4.' in filename) or ('PTV_70.4_4mm' in filename) or ('PTV_alta' in filename) or ('PTV_76' in filename) or ('PTV7' in filename) or ('PTV_high' in filename) or ('PTV_High' in filename) or ('PTV_HD.' in filename)):
            PTV_HD_path = filename
                        
    PTV_HD = sitk.ReadImage(PTV_HD_path) #read PTV HD
    
    #%%Mask Dose prescription in PTV LD
    DP_noBoneCT = generate_mask(DP_noBoneCT, PTV_LD)
    sitk.WriteImage(DP_noBoneCT, prsc_dir +'/DP_noBone_inCT.nii')
    
    #%%Generate HD and LD margins
    
    Valid_presc = DP_noBoneCT>0
    valid_in_HD = generate_mask(Valid_presc, PTV_HD)
    HD_margin = (PTV_HD-valid_in_HD)>0
    margin = (PTV_LD-PTV_HD)>0
    valid_in_LD = generate_mask(Valid_presc, margin)
    LD_margin = (margin-valid_in_LD)>0
    sitk.WriteImage(HD_margin, prsc_dir+'/HD_margin.nii')
    sitk.WriteImage(LD_margin, prsc_dir+'/LD_margin.nii')
    
    #%%Mask inverse dose prescriptions only in voxels with valid dose prescription     
    Inv_DP_noBoneCT_inCTV = generate_mask(Inv_DP_noBoneCT, Valid_presc)
    
    #%%Assign 74 Gy to all voxels in HD_margins, and 54 Gy to all voxels in LD_margins
    DP_noBoneCT_inPTV = set_mask_value(DP_noBoneCT, HD_margin, 70.4)  
    DP_noBoneCT_inPTV = set_mask_value(DP_noBoneCT_inPTV, LD_margin, 39.6)
    
    Inv_DP_noBoneCT_inPTV = set_mask_value(Inv_DP_noBoneCT_inCTV, HD_margin, (MAX_dose_noBone-70.4))
    Inv_DP_noBoneCT_inPTV = set_mask_value(Inv_DP_noBoneCT_inPTV, LD_margin, (MAX_dose_noBone-39.6))
    
    #%%Write dose prescriptions and inverse dose prescriptions        
    sitk.WriteImage(DP_noBoneCT_inPTV, prsc_dir +'/DP_noBone_CfBoost.nii')
    sitk.WriteImage(Inv_DP_noBoneCT_inPTV, prsc_dir +'/Inv_DP_noBone_CfBoost.nii')
    
    #%%Convert inverse dose prescriptions into dcm RT dose files using CT dicom series as reference dicom
    
    #First we need to flip the inverse prescriptins wrt y axis, as the dii2nii algo used originally to convert all the images flipped this wrt to original dicom images
    # Inv_DP_noBoneCT_inPTV = flip_image(Inv_DP_noBoneCT_inPTV, 'y')
    ct_ref = dcm_dir+current+'/CT/' #Path to directory of CT dicom series
    
    convert_nii_to_dicom_RTdosefile(Inv_DP_noBoneCT_inPTV, ct_ref, output_directory=prsc_dir, out_filename="Inv_DP_noBone_CfBoost.dcm")
