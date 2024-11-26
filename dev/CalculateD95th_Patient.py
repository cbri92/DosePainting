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
import matplotlib.pyplot as plt

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
Results = pd.ExcelWriter(data_supradir +'D95th_calculation.xlsx')

#Create an empty dataframe to populate as going through the loop
df = pd.DataFrame(columns=['Subject_ID', 'D95th in GTV', 'D95th in CTV', 'D99th in GTV', 'D99th in CTV','LC'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/' #Set path to subject directory
    
    #%%Read input images
    
    cell = sitk.ReadImage(subj_dir+'MRI/baseline/micro/cellApp_CORRECTED_CTV.nii') #Cellularity map
    GTV = sitk.ReadImage(subj_dir+'MRI/baseline/micro/gtv_voxel_ok_3D.nii') #GTV
    CTV = sitk.ReadImage(subj_dir+'MRI/baseline/micro/ctv_voxel_ok_3D.nii') #CTV
    
    #%%Transform cellularity map into number of cells (N0) map
    
    cell = cell*(10**8) #Convert cell map to Units: 10^4 cm-3
    orig_shape = sitk.GetArrayFromImage(cell).shape #Save the shape of the original cell map image
    voxel_volume = (cell.GetSpacing()[0]*cell.GetSpacing()[1]*cell.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
    N0 = cell*voxel_volume #Multiply the values of cellular density per voxel by the volume of one voxel to obtain the number of cells per voxel
    
    #%%Obtain vectors of N0 in voxels of GTV and CTV with valid cellularity predictions
    
    N0_GTV = allVoxInt(N0, GTV)
    N0_CTV = allVoxInt(N0, CTV)
    
    #%%Find the patient specific Dmax corresponding to the 99th and 95th percentile of the TCP
    
    #For GTV
    alpha=np.full(len(N0_GTV), a, dtype=np.float32)
    beta=np.full(len(N0_GTV), b, dtype=np.float32)
    
    # for graphing below only
    x_vals = np.arange(0,100)
    #TCP Poisson model
    vals=[np.prod( np.exp( - N0_GTV * np.exp(-alpha*x-(beta*(x**2))/n) )) for x in x_vals]
    
    TCP = np.array(vals)
    VAL_max = np.isclose(TCP, 0.99, rtol=1e-02)
    MAX_DOSE = x_vals[VAL_max][0]
    print('MAX DOSE', MAX_DOSE)
    MIN_DOSE = x_vals[TCP>0.95][0]
    print('MIN DOSE', MIN_DOSE)
    
    #Plot the TCP vs Dose
    plt.plot(x_vals,TCP)
    plt.plot([MAX_DOSE,MAX_DOSE],[0,1])
    plt.plot([MIN_DOSE,MIN_DOSE],[0,1])
    plt.legend(['TCP','99th perc dose value: '+str(MAX_DOSE), '95th perc dose value: '+str(MIN_DOSE)], loc ="upper left")
    plt.xlabel("Dose (Gy)")
    plt.ylabel("TCP")
    plt.savefig(subj_dir+'/TCP_GTV_dose_plot.jpg', bbox_inches='tight')
    plt.show()
    
    
    #For CTV       
    alpha=np.full(len(N0_CTV), a, dtype=np.float32)
    beta=np.full(len(N0_CTV), b, dtype=np.float32)
    
    # for graphing below only
    x_vals = np.arange(0,100)
    #TCP Poisson model
    vals=[np.prod( np.exp( - N0_CTV * np.exp(-alpha*x-(beta*(x**2))/n) )) for x in x_vals]
    
    TCP = np.array(vals)
    VAL_max = np.isclose(TCP, 0.99, rtol=1e-02)
    MAX_DOSE_CTV = x_vals[VAL_max][0]
    print('MAX DOSE', MAX_DOSE_CTV)
    MIN_DOSE_CTV = x_vals[TCP>0.95][0]
    print('MIN DOSE', MIN_DOSE_CTV)
    
    #Plot the TCP vs Dose
    plt.plot(x_vals,TCP)
    plt.plot([MAX_DOSE_CTV,MAX_DOSE_CTV],[0,1])
    plt.plot([MIN_DOSE_CTV,MIN_DOSE_CTV],[0,1])
    plt.legend(['TCP','99th perc dose value: '+str(MAX_DOSE_CTV), '95th perc dose value: '+str(MIN_DOSE_CTV)], loc ="upper left")
    plt.xlabel("Dose (Gy)")
    plt.ylabel("TCP")
    plt.savefig(subj_dir+'/TCP_CTV_dose_plot.jpg', bbox_inches='tight')
    plt.show()
    
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
    sheet_df = {'Subject_ID': current, 'D95th in GTV': MIN_DOSE, 'D95th in CTV':MIN_DOSE_CTV, 'D99th in GTV':MAX_DOSE, 'D99th in CTV':MAX_DOSE_CTV, 'LC':LC}
    df1 = pd.DataFrame(sheet_df, index=[0], columns=['Subject_ID', 'D95th in GTV', 'D95th in CTV', 'D99th in GTV', 'D99th in CTV', 'LC']) 
    df = df.append(df1, ignore_index=True)
    
#%%Sava data to excel spreadsheet    
df.to_excel(Results, sheet_name='Dose 95th and 99th results', index=False)
Results.save()    
    

    