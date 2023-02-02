# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 10:03:05 2023

@author: Caterina Brighi

This script register CT data to DWI data passing through T2 data in chordoma patients.
"""

#%% Import functions 

import SimpleITK as sitk
import glob
import os
import pandas as pd
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import *


#%% Set Working directory
        
data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/SacralChordoma_CNAO/nifti/' #Set working directory
# data_supradir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/test/nifti2/'

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    print('Registering images for patient: '+current)
    subj_dir = data_supradir+current #Define path to patient directory
    
    CT = sitk.ReadImage(subj_dir+'/ct.nii', sitk.sitkFloat32) #Read CT image
    #Aggiungi trova GTV, CTV, PTV  
    
    MRI_dir = subj_dir+'/MRI/'
    studies_path = [ f.path for f in os.scandir(MRI_dir) if f.is_dir() ] #Create a list of the paths to the MRI studies directories
    studies_name = [ f.name for f in os.scandir(MRI_dir) if f.is_dir() ] #Create a list of studies names
    
    for study in studies_name:
        
        print('Registering images from session: '+study)
        study_dir = MRI_dir+study
        
        if not os.listdir(study_dir) :
            print("Study "+study+" directory is empty")
        else:           
            #Create a directory to store the registration matrices
            if not os.path.exists(study_dir+'/regFiles'):#if it does not already exist, create a directory where the registration files will be saved
                os.mkdir(study_dir+'/regFiles')
            reg_dir = study_dir+'/regFiles'
            
            #Create a directory to store the registered images
            if not os.path.exists(study_dir+'/regImages'):#if it does not already exist, create a directory where the registered images will be saved
                os.mkdir(study_dir+'/regImages')
            regIm_dir = study_dir+'/regImages'
            
            #Create a directory to store the MRI info files
            if not os.path.exists(study_dir+'/info/'):#if it does not already exist, create a directory where the info files will be saved
                os.mkdir(study_dir+'/info/')
            info_dir = study_dir+'/info/'
            
            #Move all the .bval, .bvec and .json files into the info folder
            for filename in glob.glob(study_dir +'/' +'*.bval'):
                shutil.move(filename, info_dir)
                
            for filename in glob.glob(study_dir +'/' +'*.bvec'):
                shutil.move(filename, info_dir)
                
            for filename in glob.glob(study_dir +'/' +'*.json'):
                shutil.move(filename, info_dir)
                
            #%%Registration of CT image to T2 MRI image
            
            #Read T2 transversal image  
            x=glob.glob(study_dir+'/'+'*t2_tse_tra_384_*')
            if x==[]:
                print('No t2_tse_tra exists for session '+study)
            else:
                T2_path = x[0]
                T2 = sitk.ReadImage(T2_path, sitk.sitkFloat32)
                
                if os.path.exists(reg_dir+'/xf_CTonT2_1.txt'):
                    print('Applying initial transform CT-->T2')
                    tfm1 = sitk.ReadTransform(reg_dir+'/xf_CTonT2_1.txt')
                    CTinT2_temp = sitk.Resample(CT, T2, tfm1, sitk.sitkLinear) #Apply initial transform 1 to CT
                    # sitk.WriteImage(CTinT2_temp, regIm_dir+'/CTonT2_temp.nii') #Save registered image in regImages directory
                    tfm2 = sitk.ReadTransform(reg_dir+'/xf_CTonT2_2.txt')
                    CTinT2 = sitk.Resample(CTinT2_temp, T2, tfm2, sitk.sitkLinear) #Apply transform to CT
                    sitk.WriteImage(CTinT2, regIm_dir+'/CTonT2.nii') #Save registered image in regImages directory

                
                else:
                    print('Initial transform CT-->T2 does not exist')
                
                
                #%%Registration of T2 MRI image to DWI MRI image
                
                #Read DWI image
                for filename in glob.glob(study_dir+'/'+'*ep2d_diff*'):
                    if (('nii' in filename) and ('adc' not in filename)):
                        DWI_path = filename
                
                dwi = sitk.ReadImage(DWI_path, sitk.sitkFloat32) #need to only get 1 of the three images
                DWI = subsample_image(dwi,'t',1)
                
                #Resample T2 in DWI space
                T2inDWI = Resample_image(T2, DWI) #Resample in DWI space
                sitk.WriteImage(T2inDWI, regIm_dir+'/T2onDWI.nii') #Save registered image in regImages directory
                
                # #Apply rigid body registration to register T2 image to DWI image
                # T2inDWI, T2inDWI_tfm = RegisterResample_image(DWI, T2, 'same', 'Multires')
                # sitk.WriteImage(T2inDWI, regIm_dir+'/T2onDWI.nii') #Save registered image in regImages directory
                # sitk.WriteTransform(T2inDWI_tfm, reg_dir+'/xf_T2onDWI.txt') #Save registration matrix in regFiles directory
                
                #%%Resample CT in DWI
                CTinDWI = Resample_image(CTinT2, DWI) #Resample in DWI space
                sitk.WriteImage(CTinDWI, regIm_dir+'/CTonDWI.nii') #Save registered image in regImages directory
                
