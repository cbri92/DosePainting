# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 11:13:47 2023

@author: Caterina Brighi

This script:
    1. calculates the Quality Factor (QF) between the dose planned and the dose prescribed in the GTV and CTV
    2. generates QF maps in the GTV and CTV
    3. estimates the Tumour Control Probability in the GTV and CTV
"""

#%% Import functions 

import SimpleITK as sitk
import numpy as np
import os
import pandas as pd
import glob
from ImageAnalysisFunctions import flip_image, Resample_image, generate_QF_map
from ImageStatisticsFunctions import allVoxInt, calculate_QF

def calc_TCP(N0, dose, alpha, beta, n):
    
    '''This function calculates the tumour control probability given cellularity, dose, alpha, beta and n of fractions.
    N0: flattened array of cellularity values
    dose: flattened array of RBE dose values
    alpha: radiosensitivity parameter for that specific ion and cell line
    beta: radiosensitivity parameter for that specific ion and cell line
    n: number of fractions'''
    
    return np.prod( np.exp( - N0 * np.exp(-alpha*dose-(beta*(dose**2))/n)))

#Alpha and beta photons and n of fractions: we used a and b photons as per dose we are using the RBE dose
a_b = 2.4 #Gy
alpha = 0.1 #Gy^-1
beta = alpha/a_b
n = 37

#%% Set Working directory
        
data_supradir = 'path/to/patients/data/supra/directory/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

Results = {'Quality Factor': pd.DataFrame(columns=['Subject_ID','QF BL GTV', 'QF DP GTV','QF BL CTV', 'QF DP CTV']),'Tumour Control Probability': pd.DataFrame(columns=['Subject_ID','TCP BL GTV', 'TCP DP GTV','TCP BL CTV', 'TCP DP CTV']) }

ResultsWriter = pd.ExcelWriter(data_supradir +'QF_TCP_results.xlsx', engine='xlsxwriter')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current
    
    if os.path.exists(subj_dir+'/PLANS/'):
        
        rtstruct_dir = subj_dir+'/RTSTRUCT'
        prsc_dir = subj_dir+'/RTPRESC'
        
        #Read CTV HD
        for filename in glob.glob(rtstruct_dir+'/CTV*'):
            if (('CTV_7' in filename) or ('CTV_high' in filename) or ('CTV_HD' in filename)):
                CTV_HD_path = filename
                            
        CTV = sitk.ReadImage(CTV_HD_path) #read PTV LD
           
        #Read GTV
        GTV = sitk.ReadImage(rtstruct_dir+'/GTV.nii')
      
        print('Reading images for '+current) 
    
        DP_dose_path = glob.glob(subj_dir+'/PLANS/Dose_painting_4B/*RTDOSE_EFFECTIVE_PLAN*')   
        DP_dose = sitk.ReadImage(DP_dose_path) #Read DP dose
        DP_dose=DP_dose[:,:,:,0]
        DP_presc = sitk.ReadImage(prsc_dir+'/Dpainted.nii') #Read DP prescription in CTV
        DP_presc = sitk.Cast(DP_presc, sitk.sitkFloat32)
        DP_dose.SetOrigin(DP_presc.GetOrigin())
        DP_dose.SetDirection(DP_presc.GetDirection())
        DP_dose = flip_image(DP_dose, 'y')
        sitk.WriteImage(DP_dose, subj_dir+'/PLANS/Dose_painting_4B/DP_4B_reorient.nii')
        
        BL_dose_path = glob.glob(subj_dir+'/PLANS/Baseline_4B/*RTDOSE_EFFECTIVE_PLAN*')
        BL_dose = sitk.ReadImage(BL_dose_path) #Read BL dose
        BL_dose = BL_dose[:,:,:,0]
        BL_presc = sitk.ReadImage(prsc_dir+'/Boost_presc_inCTV_HD.nii') #Read BL prescription in CTV
        BL_presc = sitk.Cast(BL_presc, sitk.sitkFloat32)
        # BL_dose.SetOrigin(BL_presc.GetOrigin())
        # BL_dose.SetDirection(BL_presc.GetDirection())
        # BL_dose = flip_image(BL_dose, 'y')
        BL_dose = Resample_image(BL_dose, BL_presc, sitk.sitkLinear)
        sitk.WriteImage(BL_dose, subj_dir+'/PLANS/Baseline_4B/Baseline_4B.nii')
           
        # #Make RT struct as binary labels
        # GTV = GTV>0
        # CTV = CTV>0
    
        # #Set same origin for targets and dose grids
        # GTV.SetOrigin(DP_dose.GetOrigin())
        # CTV.SetOrigin(DP_dose.GetOrigin())
        
        #%%Generate Quality Factor Images and calculating QF
           
        print('Generating QF and QF images for '+current)   
        
        #Calculate QF for DP Plans
        QF_DP_GTV = calculate_QF(DP_presc, DP_dose, GTV)
        QF_DP_CTV = calculate_QF(DP_presc, DP_dose, CTV)
        
        #Calculate QFfor STD Plans
        QF_BL_GTV = calculate_QF(BL_presc, BL_dose, GTV)
        QF_BL_CTV = calculate_QF(BL_presc, BL_dose, CTV)
        
        #Generate QF images for DP Plans
        QF_DP_GTV_img = generate_QF_map(DP_presc, DP_dose, GTV)
        QF_DP_CTV_img = generate_QF_map(DP_presc, DP_dose, CTV)
        
        #Generate QF images for STD Plans
        QF_BL_GTV_img = generate_QF_map(BL_presc, BL_dose, GTV)
        QF_BL_CTV_img = generate_QF_map(BL_presc, BL_dose, CTV)
        
        #Make a result directory
        if not os.path.exists(subj_dir+'/Results_QF_TCP_analysis_4b'):#if it does not already exist, create a directory where the dose prescription files will be saved
            os.mkdir(subj_dir+'/Results_QF_TCP_analysis_4b')
        result_dir = subj_dir+'/Results_QF_TCP_analysis_4b/'
        
        #Save images
        
        sitk.WriteImage(QF_DP_GTV_img, result_dir+'QF_DP_GTV.nii')
        sitk.WriteImage(QF_DP_CTV_img, result_dir+'QF_DP_CTV.nii')
        
        sitk.WriteImage(QF_BL_GTV_img, result_dir+'QF_BL_GTV.nii')
        sitk.WriteImage(QF_BL_CTV_img, result_dir+'QF_BL_CTV.nii')
        
        #%%Calculate TCP within target volumes
        
        print('Calculating TCP within target volumes for '+current)
        
        reg_dir = subj_dir+'/regFiles'
        
        #Read CT image
        CT = sitk.ReadImage(subj_dir+'/ct.nii')
        
        #Read cellularity map
        cell_map = sitk.ReadImage(subj_dir+'/MRI/baseline/micro/cellApp_CORRECTED_CTV.nii') #Read cellular density map. Units: um-3
        cell_map = cell_map*(10**8) #Convert cell map to Units: 10^4 cm-3
        voxel_volume = (cell_map.GetSpacing()[0]*cell_map.GetSpacing()[1]*cell_map.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
        N0 = cell_map*voxel_volume
        
         #%%Read transformation files
        transform1 = sitk.ReadTransform(reg_dir+'/xf_CTonT2.txt')
        versor1 = sitk.VersorRigid3DTransform(transform1)
        
        transform2 = sitk.ReadTransform(reg_dir+'/xf_T2onDWI.txt')
        versor2 = sitk.VersorRigid3DTransform(transform2)
    
        composite_transform = sitk.CompositeTransform([versor2, versor1])
        inverse_transform = composite_transform.GetInverse()
        
        #%%Apply inverse transform to Dose prescriptions to bring them in CT space from DWI space
        N0 = sitk.Resample(N0, CT, inverse_transform, sitk.sitkNearestNeighbor)
    
        #%%Calculate TCP in GTV
        DP_plan_GTV = allVoxInt(DP_dose, GTV)
        BL_plan_GTV = allVoxInt(BL_dose, GTV)
        N0_GTV = allVoxInt(N0, GTV)
    
        DP_TCP_GTV = calc_TCP(N0_GTV, DP_plan_GTV, alpha, beta, n)
        BL_TCP_GTV = calc_TCP(N0_GTV, BL_plan_GTV, alpha, beta, n)
        
        #Calculate TCP in CTV
        DP_plan_CTV = allVoxInt(DP_dose, CTV)
        BL_plan_CTV = allVoxInt(BL_dose, CTV)
        N0_CTV = allVoxInt(N0, CTV)
    
        DP_TCP_CTV = calc_TCP(N0_CTV, DP_plan_CTV, alpha, beta, n)
        BL_TCP_CTV = calc_TCP(N0_CTV, BL_plan_CTV, alpha, beta, n)
        
        #%%Append QF and TCP values to Results dataframe
        
        Results['Quality Factor'] = Results['Quality Factor'].append({'Subject_ID': current,'QF BL GTV':QF_BL_GTV, 'QF DP GTV':QF_DP_GTV,'QF BL CTV':QF_BL_CTV, 'QF DP CTV':QF_DP_CTV}, ignore_index=True)    
        Results['Tumour Control Probability'] = Results['Tumour Control Probability'].append({'Subject_ID': current,'TCP BL GTV':BL_TCP_GTV, 'TCP DP GTV':DP_TCP_GTV,'TCP BL CTV':BL_TCP_CTV, 'TCP DP CTV':DP_TCP_CTV}, ignore_index=True)
        
#%%Save all dataframes to excel files here
print('Save all results to excel files')

for name, df in Results.items():
    df.to_excel(ResultsWriter, sheet_name=name, index=False)
ResultsWriter.save()
    
