# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:24:54 2023

@author: cbri3325
"""

import SimpleITK as sitk
import pandas as pd
from ImageStatisticsFunctions import allVoxInt
import numpy as np
import os

#%%Set alpha and beta particle-cell line specific values

#Set alpha/beta photons
alpha_beta_x = 2.4 #Gy for chordomas
alpha_x = 0.1 #Gy-1
beta_x = alpha_x/alpha_beta_x #Gy-2

#Values RBEmax and RBEmin for proton on chordoma cells taken from paper Paganetti, International Journal of Radiation Oncology*Biology*Physics, 2022, 112(1), 222-236.

RBEmax = 1.59
RBEmin = 1.18

a = (RBEmax*alpha_x)
b = (beta_x*(RBEmin**2))

n=37 #Number of fractions of proton treatment

#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SCB_Tutti/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('AIRC24946_R001')

#Generate excel file containing TCP results
Results = pd.ExcelWriter(data_supradir +'TCP_calculation.xlsx')

#Create an empty dataframe to populate as going through the loop
df = pd.DataFrame(columns=['Subject_ID', 'TCP in GTV', 'TCP in CTV', 'LC'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/' #Set path to subject directory
    
    #%%Read input images
    
    RBE = sitk.ReadImage(subj_dir+'RTDOSE/RBEdoseinDWI.nii') #RBE dose map
    cell = sitk.ReadImage(subj_dir+'MRI/baseline/micro/cellApp_CORRECTED_CTV.nii') #Cellularity map
    GTV = sitk.ReadImage(subj_dir+'MRI/baseline/micro/gtv_voxel_ok_3D.nii') #GTV
    CTV = sitk.ReadImage(subj_dir+'MRI/baseline/micro/ctv_voxel_ok_3D.nii') #CTV
    
    #%%Transform cellularity map into number of cells (N0) map
    
    cell = cell*(10**8) #Convert cell map to Units: 10^4 cm-3
    orig_shape = sitk.GetArrayFromImage(cell).shape #Save the shape of the original cell map image
    voxel_volume = (cell.GetSpacing()[0]*cell.GetSpacing()[1]*cell.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
    N0 = cell*voxel_volume #Multiply the values of cellular density per voxel by the volume of one voxel to obtain the number of cells per voxel
    
    #%%Obtain vectors of RBE and N0 in voxels of GTV and CTV with valid cellularity predictions
    
    RBE_GTV = allVoxInt(RBE, GTV)
    RBE_CTV = allVoxInt(RBE, CTV)
    
    N0_GTV = allVoxInt(N0, GTV)
    N0_CTV = allVoxInt(N0, CTV)
    
    #%%Calculate TCP in GTV
       
    alpha=np.full(len(N0_GTV), a, dtype=np.float32)
    beta=np.full(len(N0_GTV), b, dtype=np.float32)

    TCP_GTV = np.prod( np.exp( - N0_GTV * np.exp(-alpha*RBE_GTV-(beta*(RBE_GTV**2))/n)))
    
    #%%Calculate TCP in CTV
       
    alpha=np.full(len(N0_CTV), a, dtype=np.float32)
    beta=np.full(len(N0_CTV), b, dtype=np.float32)
    
    TCP_CTV = np.prod( np.exp( - N0_CTV * np.exp(-alpha*RBE_CTV-(beta*(RBE_CTV**2))/n)))
    
    #%%determine LC
    if current == 'AIRC24946_R012':
        LC = 0
    elif current == 'AIRC24946_R029':
        LC = 0
    elif current == 'AIRC24946_R035':
        LC = 0
    elif current == 'AIRC24946_R041':
        LC = 0
    elif current == 'AIRC24946_R048':
        LC = 0
    elif current == 'AIRC24946_R052':
        LC = 0
    else:
        LC = 1
    
#%%Save TCP results           
    sheet_df = {'Subject_ID': current, 'TCP in GTV':TCP_GTV, 'TCP in CTV':TCP_CTV, 'LC':LC}
    df1 = pd.DataFrame(sheet_df, index=[0], columns=['Subject_ID', 'TCP in GTV', 'TCP in CTV', 'LC']) 
    df = df.append(df1, ignore_index=True)
    
#%%Sava data to excel spreadsheet    
df.to_excel(Results, sheet_name='TCP results', index=False)
Results.save()    
    

    