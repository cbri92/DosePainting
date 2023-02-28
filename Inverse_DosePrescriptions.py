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
# from ConvertNiftiImage_ToDoseFiles import * #Use this functions only when using the dicom RTDOSE file as reference dicom
from ConvertNii_ToDoseFiles import * #Use this functions only when using the dicom CT series file as reference dicom

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('AIRC24946_R052')
# subjs_name=['AIRC24946_R029']
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
    # DOSE = sitk.ReadImage(subj_dir+'/RTDOSE/PHYS_v0.nii')
    
    #Read dose prescription images
    DP_origDWI = sitk.ReadImage(dose_dir +'/Dose_optimised_CTVorig_final.nii')
    DP_noBoneDWI = sitk.ReadImage(dose_dir +'/Dose_optimised_CTVnoBone_final.nii')
    
    #Read transformation files
    transform1 = sitk.ReadTransform(reg_dir+'/xf_CTonT2.txt')
    versor1 = sitk.VersorRigid3DTransform(transform1)
    
    transform2 = sitk.ReadTransform(reg_dir+'/xf_T2onDWI.txt')
    versor2 = sitk.VersorRigid3DTransform(transform2)

    composite_transform = sitk.CompositeTransform([versor2, versor1])
    inverse_transform = composite_transform.GetInverse()
    
    #Apply inverse transform to Dose prescriptions to bring them in CT space from DWI space
    DP_origCT = sitk.Resample(DP_origDWI, CT, inverse_transform, sitk.sitkNearestNeighbor)
    DP_noBoneCT = sitk.Resample(DP_noBoneDWI, CT, inverse_transform, sitk.sitkNearestNeighbor)
   
    #Calculate Max and Min dose prescription values
    Stats_dose_orig = getNonZeroStats(DP_origCT)
    MAX_dose_orig = int(Stats_dose_orig['Max intensity'])
    MIN_dose_orig = int(Stats_dose_orig['Min intensity'])
    MEAN_dose_orig = int(Stats_dose_orig['Mean intensity'])
    
    Stats_dose_noBone = getNonZeroStats(DP_noBoneCT)
    MAX_dose_noBone = int(Stats_dose_noBone['Max intensity'])
    MIN_dose_noBone = int(Stats_dose_noBone['Min intensity'])
    MEAN_dose_noBone = int(Stats_dose_noBone['Mean intensity'])
    
    #Generate inverse dose prescription as Dinv = Dmax - Dpresc, where Dmax = maximum dose estimated for each patient
    Inv_DP_origCT = MAX_dose_orig-DP_origCT
    Inv_DP_noBoneCT = MAX_dose_noBone-DP_noBoneCT
    
    #Read CTV and PTV LD
    CTV = DP_origCT>0
    
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'PTV*'):
        if (('PTV_5' in filename) or ('PTV5' in filename) or ('PTV_3' in filename) or ('PTV_low' in filename) or ('PTV_LD' in filename)):
            PTV_path = filename
        # elif 'CTV_AIRC' in filename:
        #     CTV_path = filename
                        
    # CTV = sitk.ReadImage(CTV_path) #read CTV 
    PTV = sitk.ReadImage(PTV_path) #read PTV
    margin = (PTV-CTV)>0
    # sitk.WriteImage(margin, subj_dir+'/margin.nii')
    
    #Generate mask of dose prescriptions and inverse dose prescriptions in CTV    
    DP_origCT_inCTV = generate_mask(DP_origCT, CTV)
    DP_noBoneCT_inCTV = generate_mask(DP_noBoneCT, CTV)    
    
    Inv_DP_origCT_inCTV = generate_mask(Inv_DP_origCT, CTV)
    Inv_DP_noBoneCT_inCTV = generate_mask(Inv_DP_noBoneCT, CTV)
    
    #Assign min dose prescription to all voxels in margins CTV-PTV, where min dose is the minimum dose estimated for each patient
    DP_origCT_inPTV = set_mask_value(DP_origCT_inCTV, margin, MEAN_dose_orig)
    DP_noBoneCT_inPTV = set_mask_value(DP_noBoneCT_inCTV, margin, MEAN_dose_noBone)    
    
    Inv_DP_origCT_inPTV = set_mask_value(Inv_DP_origCT_inCTV, margin, (MAX_dose_orig-MEAN_dose_orig))
    Inv_DP_noBoneCT_inPTV = set_mask_value(Inv_DP_noBoneCT_inCTV, margin, (MAX_dose_noBone-MEAN_dose_noBone))
    
    # #Resample dose prescriptions in dose space
    # DP_origDOSE_inPTV = Resample_image(DP_origCT_inPTV, DOSE, sitk.sitkNearestNeighbor)
    # DP_noBoneDOSE_inPTV = Resample_image(DP_noBoneCT_inPTV, DOSE, sitk.sitkNearestNeighbor)
    
    # Inv_DP_origDOSE_inPTV = Resample_image(Inv_DP_origCT_inPTV, DOSE, sitk.sitkNearestNeighbor)
    # Inv_DP_noBoneDOSE_inPTV = Resample_image(Inv_DP_noBoneCT_inPTV, DOSE, sitk.sitkNearestNeighbor)
    
    #Write dose prescriptions and inverse dose prescriptions        
    sitk.WriteImage(DP_origCT_inPTV, prsc_dir +'/DP_orig.nii')
    sitk.WriteImage(DP_noBoneCT_inPTV, prsc_dir +'/DP_noBone.nii')
    
    sitk.WriteImage(Inv_DP_origCT_inPTV, prsc_dir +'/Inv_DP_orig.nii')
    sitk.WriteImage(Inv_DP_noBoneCT_inPTV, prsc_dir +'/Inv_DP_noBone.nii')
    
    # sitk.WriteImage(DP_origDOSE_inPTV, prsc_dir +'/DP_orig.nii')
    # sitk.WriteImage(DP_noBoneDOSE_inPTV, prsc_dir +'/DP_noBone.nii')
    
    # sitk.WriteImage(Inv_DP_origDOSE_inPTV, prsc_dir +'/Inv_DP_orig.nii')
    # sitk.WriteImage(Inv_DP_noBoneDOSE_inPTV, prsc_dir +'/Inv_DP_noBone.nii')
    
    #Convert inverse dose prescriptions into dcm RT dose files using CT dicom series as reference dicom
        #First we need to flip the inverse prescriptins wrt y axis, as the dii2nii algo used originally to convert all the images flipped this wrt to original dicom images
    Inv_DP_origCT_inPTV = flip_image(Inv_DP_origCT_inPTV, 'y')
    Inv_DP_noBoneCT_inPTV = flip_image(Inv_DP_noBoneCT_inPTV, 'y')

    ct_ref = subj_dir+'/CT/' #Path to directory of CT dicom series
    convert_nii_to_dicom_RTdosefile(Inv_DP_origCT_inPTV, ct_ref, output_directory=prsc_dir, out_filename="Inv_DP_orig.dcm")
    convert_nii_to_dicom_RTdosefile(Inv_DP_noBoneCT_inPTV, ct_ref, output_directory=prsc_dir, out_filename="Inv_DP_noBone.dcm")

    # #Convert inverse dose prescriptions into dcm RT dose files using RTDOSE dicom file as reference dicom
    # dcm_ref = glob.glob(subj_dir+'/RTDOSE_dcm_ref/'+'*.dcm')[0]
    # convert_nifti_to_dicom_RTdosefile(Inv_DP_origDOSE_inPTV, dcm_ref, output_directory=prsc_dir, out_filename="Inv_DP_orig.dcm")
    # convert_nifti_to_dicom_RTdosefile(Inv_DP_noBoneDOSE_inPTV, dcm_ref, output_directory=prsc_dir, out_filename="Inv_DP_noBone.dcm")
