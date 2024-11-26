# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 14:22:51 2023

@author: cbri3325
"""


#%% Import functions 

import SimpleITK as sitk
import os
import glob
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import allVoxInt

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SkullBaseChordoma_CNAO/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names


#Generate excel file containing voxelwise results
Results = pd.ExcelWriter(data_supradir +'Recurrence_cellularity_analysis.xlsx')

#Create an empty dataframe to populate as going through the loop
df = pd.DataFrame(columns=['Subject_ID', '% Vol of recurrence in GTV', 'Mean cellularity in recurrent GTV', 'Mean cellularity in non targeted GTV', 'Mean dose in recurrent GTV', 'Mean dose in non targeted GTV'])


#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    print('Converting cellularity map into CT space for '+current)
    subj_dir = data_supradir+current
    subj_name = current
    
    #Generate output directory where to save figures
    if not os.path.exists(subj_dir+'/Recurrence GTV analysis'):#if it does not already exist, create a directory where the figures will be saved
        os.mkdir(subj_dir+'/Recurrence GTV analysis')
    out_dir = subj_dir+'/Recurrence GTV analysis'
    
    #Read directories' path
    micro_dir = subj_dir+'/MRI/baseline/micro'
    reg_dir = subj_dir+'/regFiles'
    str_dir = subj_dir +'/RTSTRUCT'

    #Read cellularity map
    cell_map = sitk.ReadImage(micro_dir+'/cellApp_CORRECTED.nii.gz')
    
    #Read RBE dose map
    RBE = sitk.ReadImage(subj_dir+'/RTDOSE/RBEdoseinCT.nii')
    
    #Read CT image
    CT = sitk.ReadImage(subj_dir+'/ct.nii')
    
    #%%Read transformation files
    transform1 = sitk.ReadTransform(reg_dir+'/xf_CTonT2.txt')
    versor1 = sitk.VersorRigid3DTransform(transform1)
    
    transform2 = sitk.ReadTransform(reg_dir+'/xf_T2onDWI.txt')
    versor2 = sitk.VersorRigid3DTransform(transform2)

    composite_transform = sitk.CompositeTransform([versor2, versor1])
    inverse_transform = composite_transform.GetInverse()
    
    #%%Apply inverse transform to Dose prescriptions to bring them in CT space from DWI space
    cell_inCT = sitk.Resample(cell_map, CT, inverse_transform, sitk.sitkNearestNeighbor)
    sitk.WriteImage(cell_inCT, subj_dir+'/cell_density_inCT.nii')
    cell_inCT = sitk.ReadImage(subj_dir+'/cell_density_inCT.nii')
    #%%Read GTV and recurrence contours
    
    GTV = sitk.ReadImage(str_dir+'/GTV.nii')
    Rec = sitk.ReadImage(str_dir+'/Recurrence_fixOrigDir.nii')
    
    noRec_GTV = (GTV-Rec)>0
    # cell_ok = cell_inCT>0
    
    #%%Calculate % vol rec overlapping withn GTV, mean cellularity in area of overlap and in rest of GTV
    
    #Calculate % Volume of recurrence in GTV
    Tar_rec = (GTV+Rec)>1 #Volume of recurrence overlapping with GTV
    max_val = getMaxVox(Tar_rec)[0]
    if max_val<1:
        print('Recurrence does not overlay on GTV for '+current)
        
        Perc_targeted_rec = 0
        cell_trec_mean = 'NaN'
        cell_notrec_mean = 'NaN'
    
    else:
        sitk.WriteImage(Tar_rec, str_dir+'/Rec_inGTV.nii')
        
        Vol_targeted_rec = round(getVolumeRoi(Tar_rec, cell_inCT))
        Vol_Rec = round(getVolumeRoi(Rec, cell_inCT))
        Perc_targeted_rec = (Vol_targeted_rec/Vol_Rec)*100
    
        #Calculate mean cellularity in area of overlap
        cell_trec_mean = getMeanRoi(Tar_rec, cell_inCT)
        
        #Plot distribution of cellularity in targeted recurrence
        cell_trec = allVoxInt(cell_inCT, Tar_rec)
        plt.hist(cell_trec, bins=20)
        plt.xlabel("Cellular density in recurrent GTV")
        plt.ylabel("Counts")
        plt.savefig(out_dir+'/cellularity_inRecGTV.jpg')
        plt.close()
        
        #Calculate mean cellularity in rest of GTV
        cell_notrec_mean = getMeanRoi(noRec_GTV, cell_inCT)
        
        #Plot distribution of cellularity in NON recurrent GTV
        cell_trec = allVoxInt(cell_inCT, noRec_GTV)
        plt.hist(cell_trec, bins=20)
        plt.xlabel("Cellular density in NOT recurrent GTV")
        plt.ylabel("Counts")
        plt.savefig(out_dir+'/cellularity_inNORecGTV.jpg')
        plt.close()
        
        #Calculate mean dose in area of overlap
        dose_trec_mean = getMeanRoi(Tar_rec, RBE)
        
        #Calculate mean dose in rest of GTV
        dose_notrec_mean = getMeanRoi(noRec_GTV, RBE)
        
    #%%Save statistics of various VOIs            
    sheet_df = {'Subject_ID': subj_name, '% Vol of recurrence in GTV':Perc_targeted_rec, 'Mean cellularity in recurrent GTV':cell_trec_mean, 'Mean cellularity in non targeted GTV':cell_notrec_mean, 'Mean dose in recurrent GTV':dose_trec_mean, 'Mean dose in non targeted GTV':dose_notrec_mean}
    df1 = pd.DataFrame(sheet_df, index=[0], columns=['Subject_ID', '% Vol of recurrence in GTV', 'Mean cellularity in recurrent GTV', 'Mean cellularity in non targeted GTV', 'Mean dose in recurrent GTV', 'Mean dose in non targeted GTV']) 
    df = df.append(df1, ignore_index=True)
    
#%%Sava data to excel spreadsheet    
df.to_excel(Results, sheet_name='Sheet1', index=False)
Results.save()