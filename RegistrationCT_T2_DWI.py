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
    
    for filename in glob.glob(subj_dir+'/RTSTRUCT/'+'*nii.gz'):
        if (('gtv.nii.gz' in filename) or ('GTV.nii.gz' in filename) or ('GTVmrf' in filename) or ('GTV_POST_op' in filename) or ('GTV_CT' in filename) or ('GTV_ec_TAC' in filename)):
            GTV_path = filename
            print(GTV_path)
        elif (('CTV_39.6.nii.gz' in filename) or ('CTV39.6' in filename) or ('CTV_41.4' in filename) or ('CTV_43.2' in filename) or ('CTV_bassa' in filename) or ('CTV_LD.nii.gz' in filename) or ('CTVLD' in filename) or ('CTV_low' in filename)):
            CTV_path = filename
            print(CTV_path)
        elif (('PTV_39.6' in filename) or ('PTV39.6' in filename) or ('PTV_41.4' in filename) or ('PTV_43.2' in filename) or ('PTV_bassa' in filename) or ('PTV_LD.nii.gz' in filename) or ('PTVLD' in filename) or ('PTV_low' in filename)):
            PTV_path = filename
            print(PTV_path)
    
    GTV = sitk.ReadImage(GTV_path)
    CTV = sitk.ReadImage(CTV_path)  
    PTV = sitk.ReadImage(PTV_path)  
    
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
                    tfm2 = sitk.ReadTransform(reg_dir+'/xf_CTonT2_2.txt')
                    CTinT2_temp = sitk.Resample(CT, T2, tfm1, sitk.sitkLinear) #Apply initial transform 1 to CT
                    GTVinT2_temp = sitk.Resample(GTV, T2, tfm1, sitk.sitkNearestNeighbor) #Apply initial transform 1 to GTV
                    CTVinT2_temp = sitk.Resample(CTV, T2, tfm1, sitk.sitkNearestNeighbor) #Apply initial transform 1 to CTV
                    PTVinT2_temp = sitk.Resample(PTV, T2, tfm1, sitk.sitkNearestNeighbor) #Apply initial transform 1 to PTV
                    # sitk.WriteImage(CTinT2_temp, regIm_dir+'/CTonT2_temp.nii') #Save registered image in regImages directory  
                    CTinT2 = sitk.Resample(CTinT2_temp, T2, tfm2, sitk.sitkLinear) #Apply transform to CT
                    GTVinT2 = sitk.Resample(GTVinT2_temp, T2, tfm2, sitk.sitkNearestNeighbor) #Apply transform to GTV
                    CTVinT2 = sitk.Resample(CTVinT2_temp, T2, tfm2, sitk.sitkNearestNeighbor) #Apply transform to CTV
                    PTVinT2 = sitk.Resample(PTVinT2_temp, T2, tfm2, sitk.sitkNearestNeighbor) #Apply transform to PTV
                    
                    sitk.WriteImage(CTinT2, regIm_dir+'/CTonT2.nii') #Save registered image in regImages directory
                    sitk.WriteImage(GTVinT2, regIm_dir+'/GTVonT2.nii') #Save registered image in regImages directory
                    sitk.WriteImage(CTVinT2, regIm_dir+'/CTVonT2.nii') #Save registered image in regImages directory
                    sitk.WriteImage(PTVinT2, regIm_dir+'/PTVonT2.nii') #Save registered image in regImages directory
  
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
                T2inDWI = Resample_image(T2, DWI, sitk.sitkLinear) #Resample in DWI space
                sitk.WriteImage(T2inDWI, regIm_dir+'/T2onDWI.nii') #Save registered image in regImages directory
                
                # #Apply rigid body registration to register T2 image to DWI image
                # T2inDWI, T2inDWI_tfm = RegisterResample_image(DWI, T2, 'same', 'Multires')
                # sitk.WriteImage(T2inDWI, regIm_dir+'/T2onDWI.nii') #Save registered image in regImages directory
                # sitk.WriteTransform(T2inDWI_tfm, reg_dir+'/xf_T2onDWI.txt') #Save registration matrix in regFiles directory
                
                #%%Resample CT, GTV, CTV and PTV in DWI
                CTinDWI = Resample_image(CTinT2, DWI, sitk.sitkLinear) #Resample in DWI space
                GTVinDWI = Resample_roi(GTVinT2, DWI) #Resample in DWI space
                CTVinDWI = Resample_roi(CTVinT2, DWI) #Resample in DWI space
                PTVinDWI = Resample_roi(PTVinT2, DWI) #Resample in DWI space
                
                sitk.WriteImage(CTinDWI, regIm_dir+'/CTonDWI.nii') #Save registered image in regImages directory
                sitk.WriteImage(GTVinDWI, regIm_dir+'/GTVonDWI.nii') #Save registered image in regImages directory
                sitk.WriteImage(CTVinDWI, regIm_dir+'/CTVonDWI.nii') #Save registered image in regImages directory
                sitk.WriteImage(PTVinDWI, regIm_dir+'/PTVonDWI.nii') #Save registered image in regImages directory
                
                #%%Create 5 mm and 8 mm expansions to GTV in DWI space
                Exp_5mm = sitk.BinaryDilate(GTVinDWI, (round(5/DWI.GetSpacing()[0]),round(5/DWI.GetSpacing()[1]),round(5/DWI.GetSpacing()[2])))
                Exp_8mm = sitk.BinaryDilate(GTVinDWI, (round(8/DWI.GetSpacing()[0]),round(8/DWI.GetSpacing()[1]),round(8/DWI.GetSpacing()[2])))

                sitk.WriteImage(Exp_5mm, regIm_dir+'/GTV_exp5mm.nii') #Save registered image in regImages directory
                sitk.WriteImage(Exp_8mm, regIm_dir+'/GTV_exp8mm.nii') #Save registered image in regImages directory
