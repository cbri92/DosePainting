# -*- coding: utf-8 -*-
"""
Created on Tue May  2 18:19:20 2023

@author: Caterina Brighi

This script finds alphaPOP through an optimization problem, 
which consisted in finding the value of α (αPOP) for Np patients so that the population TCP (TCPPOP), 
modelled through a Poissonian statistics, matched the observed LCPOP. 

Following method published in supplementary material Step 3 of: Buizza et al. Radiotherapy and Oncology 137 (2019) 32–37.
"""

import numpy as np
import pandas as pd
import os
from scipy.optimize import minimize,Bounds, minimize_scalar, brentq
import matplotlib.pyplot as plt
from ImageStatisticsFunctions import allVoxInt
from ImageAnalysisFunctions import generate_mask
import SimpleITK as sitk

#%%Define function for optimization problem

def function(alphaPAT, N0, Dose_inTarget, tcp):   
    prod = np.prod( np.exp( - N0 * np.exp(-alphaPAT*Dose_inTarget)))
    f = tcp-prod
    return f

def Func(x):
    prod = np.prod( np.exp( - N0 * np.exp(-x*Dose_inTarget)))
    f = tcp-prod
    return f

#%%Write new excel spreadsheet

Results = pd.ExcelWriter('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/alpha_patients.xlsx')
res_df = pd.DataFrame(columns=['Patient_ID', 'alphaPAT'])

#%%Import DVH info dataframe
dvh_df = pd.read_excel('C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/DVH_info_GTV.xlsx', sheet_name='DVH stats', index_col='Patient_ID')

#%%Calculate alphaPAT for recurrent patients

data_dir='C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti/'

subjs_path = [f.path for f in os.scandir(data_dir) if f.is_dir()] #Create a list of the paths to the subjects directories
subjs_name = [f.name for f in os.scandir(data_dir) if f.is_dir()] #Create a list of subjects names
subjs_name.remove('P29')
subjs_name.remove('P38')
subjs_name.remove('P39')
subjs_name.remove('P41')
subjs_name.remove('P45')
subjs_name.remove('P50')

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name: 
    
    subj_dir = data_dir+current    
    subj_name = current
     
    print(current)

    #Read CT image
    ct = sitk.ReadImage(subj_dir + '/ct.nii')
    
    #Read cellularity image
    # cell = sitk.ReadImage(subj_dir+'/cell_density_inCT.nii')
    
    #Read RBE in CT
    RBE = sitk.ReadImage(subj_dir + '/RTDOSE/RBE_inCT.nii')
    
    # Read GTV image  
    GTV = sitk.ReadImage(subj_dir + '/GTV_ok_inCT.nii')
    GTV.SetOrigin(ct.GetOrigin())
    GTV.SetDirection(ct.GetDirection())
    
    # Calculate Dose_inTarget
    Dose_inTarget = allVoxInt(RBE, GTV)
    
    # Calculate N0
    # cell_map = cell*(10**8) #Convert cell map to Units: 10^4 cm-3
    # cell_map_GTV = generate_mask(cell_map, GTV)
    # voxel_volume = (cell_map_GTV.GetSpacing()[0]*cell_map_GTV.GetSpacing()[1]*cell_map_GTV.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
    # cell_dens = allVoxInt(cell_map_GTV, GTV) #Define array of values with cellular density>0
    # N0 = cell_dens*voxel_volume #Multiply the values of cellular density per voxel by the volume of one voxel to obtain the number of cells per voxel
    N0=10**7
    
    # Obtain TCP from spreadsheet
    tcp = dvh_df.loc[subj_name]['LC']

    #%% Find alphaPAT through optimization function
    
    x = brentq(Func, 0.05, 1.0)
    alphaPAT=x
    
    alphaPAT0_vals = np.linspace(0,2,100)
    opt = [function(alpha, N0, Dose_inTarget, tcp) for alpha in alphaPAT0_vals]
    f = np.array(opt)
    alphaPATmax = alphaPAT0_vals[f<-0.01][0]
    
    plt.plot(alphaPAT0_vals, opt)
    plt.plot([alphaPATmax,alphaPATmax],[min(opt),max(opt)])
    plt.legend(['f','alphaPAT: '+str(round(alphaPATmax,2))], loc ="lower right")
    plt.xlabel("alpha (Gy-1)")
    plt.ylabel("f")
    # plt.savefig(subj_dir+'/alphaPAT.jpg')
    plt.show()

    #Save stats in dataframe
    print('Patient '+subj_name+' , alphaPAT = '+str(alphaPATmax))
    res_dict={'Patient_ID':subj_name, 'alphaPAT':alphaPATmax}
    res_df = res_df.append(res_dict, ignore_index=True)
    
res_df.to_excel(Results, sheet_name='alphaPAT', index=False)
Results.save()