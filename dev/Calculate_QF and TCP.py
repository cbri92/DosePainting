# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 11:13:47 2023

@author: cbri3325
"""

#%% Import functions 

import SimpleITK as sitk
import numpy as np
import os
import pandas as pd
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import *

def calc_TCP(N0, dose, alpha, beta, n):
    return np.prod( np.exp( - N0 * np.exp(-alpha*dose-(beta*(dose**2))/n)))

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
# subjs_name=['AIRC24946_R012']

Results = {'Quality Factor': pd.DataFrame(columns=['Subject_ID','QF STD GTV', 'QF DP GTV','QF STD CTV', 'QF DP CTV','QF STD PTV', 'QF DP PTV']),'Tumour Control Probability': pd.DataFrame(columns=['Subject_ID','TCP STD GTV', 'TCP DP GTV','TCP STD CTV', 'TCP DP CTV','TCP STD PTV', 'TCP DP PTV']) }

ResultsWriter = pd.ExcelWriter(data_supradir +'QF_TCP_results.xlsx', engine='xlsxwriter')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current
    
    rtstruct_dir = subj_dir+'/RTSTRUCT'
    prsc_dir = subj_dir+'/RTPRESC'
    
    #Read PTV LD
    for filename in glob.glob(rtstruct_dir+'/PTV*'):
        if (('PTV_5' in filename) or ('PTV5' in filename) or ('PTV_3' in filename) or ('PTV_low' in filename) or ('PTV_LD' in filename)):
            PTV_LD_path = filename
                        
    PTV = sitk.ReadImage(PTV_LD_path) #read PTV LD
    
    #Read CTV LD
    for filename in glob.glob(rtstruct_dir+'/CTV*'):
        if (('CTV_5' in filename) or ('CTV5' in filename) or ('CTV_3' in filename) or ('CTV_low' in filename) or ('CTV_LD' in filename)):
            CTV_LD_path = filename
                        
    CTV = sitk.ReadImage(CTV_LD_path) #read CTV LD
    
    #Read GTV
    GTV = sitk.ReadImage(rtstruct_dir+'/GTV.nii')
  
    print('Reading images for '+current) 

    DP_dose_path = glob.glob(subj_dir+'/Dose_painted/*RTDOSE_EFFECTIVE_PLAN*')

    DP_dose = sitk.ReadImage(DP_dose_path) #Read DP dose in PTV
    DP_dose=DP_dose[:,:,:,0]
    DP_presc = sitk.ReadImage(prsc_dir+'/DP_noBone_CfBoost.nii') #Read DP prescription in PTV
    DP_presc = sitk.Cast(DP_presc, sitk.sitkFloat32)
    DP_dose.SetOrigin(DP_presc.GetOrigin())
    DP_dose.SetDirection(DP_presc.GetDirection())
    DP_dose = flip_image(DP_dose, 'y')
    sitk.WriteImage(DP_dose, subj_dir+'/Dose_painted/DP_dose_reorient.nii')
    # sitk.WriteImage(Rel_dose, subj_dir+'/Dose_painted/Rel_dose.nii')
    
    STD_dose = sitk.ReadImage(subj_dir+'/RTDOSE/RBEdoseinCT.nii') #Read STD dose
    STD_presc = sitk.ReadImage(prsc_dir+'/Boost_presc_inPTV.nii') #Read STD prescription in PTV
    STD_presc = sitk.Cast(STD_presc, sitk.sitkFloat32)
       
    # #Make RT struct as binary labels
    # GTV = GTV>0
    # CTV = CTV>0
    # PTV = PTV>0

    # #Set same origin for targets and dose grids
    # GTV.SetOrigin(DP_dose.GetOrigin())
    # CTV.SetOrigin(DP_dose.GetOrigin())
    # PTV.SetOrigin(DP_dose.GetOrigin())
    
    #%%Generate Quality Factor Images and calculating QF
       
    print('Generating QF and QF images for '+current)   
    
    #Calculate QF for DP Plans
    QF_DP_GTV = calculate_QF(DP_presc, DP_dose, GTV)
    QF_DP_CTV = calculate_QF(DP_presc, DP_dose, CTV)
    QF_DP_PTV = calculate_QF(DP_presc, DP_dose, PTV)
    
    #Calculate QFfor STD Plans
    QF_STD_GTV = calculate_QF(STD_presc, STD_dose, GTV)
    QF_STD_CTV = calculate_QF(STD_presc, STD_dose, CTV)
    QF_STD_PTV = calculate_QF(STD_presc, STD_dose, PTV)
    
    #Generate QF images for DP Plans
    QF_DP_GTV_img = generate_QF_map(DP_presc, DP_dose, GTV)
    QF_DP_CTV_img = generate_QF_map(DP_presc, DP_dose, CTV)
    QF_DP_PTV_img = generate_QF_map(DP_presc, DP_dose, PTV)
    
    #Generate QF images for STD Plans
    QF_STD_GTV_img = generate_QF_map(STD_presc, STD_dose, GTV)
    QF_STD_CTV_img = generate_QF_map(STD_presc, STD_dose, CTV)
    QF_STD_PTV_img = generate_QF_map(STD_presc, STD_dose, PTV)
    
    #Make a result directory
    if not os.path.exists(subj_dir+'/Results_QF_TCP_analysis'):#if it does not already exist, create a directory where the dose prescription files will be saved
        os.mkdir(subj_dir+'/Results_QF_TCP_analysis')
    result_dir = subj_dir+'/Results_QF_TCP_analysis/'
    
    #Save images
    
    sitk.WriteImage(QF_DP_GTV_img, result_dir+'QF_DP_GTV.nii')
    sitk.WriteImage(QF_DP_CTV_img, result_dir+'QF_DP_CTV.nii')
    sitk.WriteImage(QF_DP_PTV_img, result_dir+'QF_DP_PTV.nii')
    
    sitk.WriteImage(QF_STD_GTV_img, result_dir+'QF_Boost_GTV.nii')
    sitk.WriteImage(QF_STD_CTV_img, result_dir+'QF_Boost_CTV.nii')
    sitk.WriteImage(QF_STD_PTV_img, result_dir+'QF_Boost_PTV.nii')

    
    #%%Calculate TCP within target volumes
    
    print('Calculating TCP within target volumes for '+current)
    
    reg_dir = subj_dir+'/regFiles'
    
    #Read CT image
    CT = sitk.ReadImage(subj_dir+'/ct.nii')
    
    #Read cellularity map
    cell_map = sitk.ReadImage(subj_dir+'/MRI/baseline/orig/cellApp_CORRECTED.nii') #Read cellular density map. Units: um-3
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
    STD_plan_GTV = allVoxInt(STD_dose, GTV)
    N0_GTV = allVoxInt(N0, GTV)

    DP_TCP_GTV = calc_TCP(N0_GTV, DP_plan_GTV, 0.159, 0.058, 37)
    STD_TCP_GTV = calc_TCP(N0_GTV, STD_plan_GTV, 0.159, 0.058, 37)
    
    #Calculate TCP in PTV
    DP_plan_CTV = allVoxInt(DP_dose, CTV)
    STD_plan_CTV = allVoxInt(STD_dose, CTV)
    N0_CTV = allVoxInt(N0, CTV)

    DP_TCP_CTV = calc_TCP(N0_CTV, DP_plan_CTV, 0.159, 0.058, 37)
    STD_TCP_CTV = calc_TCP(N0_CTV, STD_plan_CTV, 0.159, 0.058, 37)
    
    #Calculate TCP in PTV
    DP_plan_PTV = allVoxInt(DP_dose, PTV)
    STD_plan_PTV = allVoxInt(STD_dose, PTV)
    N0_PTV = allVoxInt(N0, PTV)

    DP_TCP_PTV = calc_TCP(N0_PTV, DP_plan_PTV, 0.159, 0.058, 37)
    STD_TCP_PTV = calc_TCP(N0_PTV, STD_plan_PTV, 0.159, 0.058, 37)
    
    #%%Append QF and TCP values to Results dataframe
    
    Results['Quality Factor'] = Results['Quality Factor'].append({'Subject_ID': current,'QF STD GTV':QF_STD_GTV, 'QF DP GTV':QF_DP_GTV,'QF STD CTV':QF_STD_CTV, 'QF DP CTV':QF_DP_CTV,'QF STD PTV':QF_STD_PTV, 'QF DP PTV':QF_DP_PTV}, ignore_index=True)    
    Results['Tumour Control Probability'] = Results['Tumour Control Probability'].append({'Subject_ID': current,'TCP STD GTV':STD_TCP_GTV, 'TCP DP GTV':DP_TCP_GTV,'TCP STD CTV':STD_TCP_CTV, 'TCP DP CTV':DP_TCP_CTV,'TCP STD PTV':STD_TCP_PTV, 'TCP DP PTV':DP_TCP_PTV}, ignore_index=True)
    
#%%Save all dataframes to excel files here
print('Save all results to excel files')

for name, df in Results.items():
    df.to_excel(ResultsWriter, sheet_name=name, index=False)
ResultsWriter.save()
    