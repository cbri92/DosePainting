# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:31:51 2022

@author: Caterina Brighi

This script is used to take information about cellularity, alpha and beta and use them
to derive patient-specific dose prescriptions by means of a TCP model (poissonian, LQ).
"""

#%% Import functions 

import SimpleITK as sitk
import os
import numpy as np

def linear_Dose_N0_mapping(N0, Dmin, Dmax):
    """This function linearly maps values of cellularity to values of dose.
        N0: array of values of cellularity
        Dmin: minimum dose to deliver to voxels
        Dmax: maximum dose to deliver to voxels"""
        
    N0_max = max(N0)
    N0_min = min(N0)        
    Dose = Dmin + ((Dmax-Dmin)*(N0-N0_min))/(N0_max-N0_min)
    return Dose

def generate_mask(image, roi):
    
    '''Returns the masked image of the roi applied to the image.
    Remember to save the masked image file after applying this function.'''
    
    masking = sitk.MaskImageFilter() 
    mask = masking.Execute(image, roi)
    return mask

def allVoxInt(image, roi):
    
    '''Returns a flatten array (in 2D) of the intensity values from all pixels in roi applied to an image.'''
    
    image = sitk.GetArrayFromImage(image)
    roi = sitk.GetArrayFromImage(roi)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    masked = image[~np.array(mask)]
    return masked
    
def set_mask_value(image, mask, value):
    
    '''Function that takes as input an image, a mask and a value, and returns the original 
    image with the values of the pixels corresponding
    to the binary mask set to the specified value.'''
    
    msk32 = sitk.Cast(mask, sitk.sitkFloat32)
    return sitk.Cast(sitk.Cast(image, sitk.sitkFloat32) *
                     sitk.InvertIntensity(msk32, maximum=1.0) + 
                     msk32*value, image.GetPixelID())

#%% Set Working directory
        
data_supradir = 'path/to/patients/data/supra/directory/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current+'/MRI/baseline/micro/' #Set path to subject directory where cell density map image and GTV image are stored
    if not os.path.exists(data_supradir+current+'/MRI/baseline/Prescribed_dose'):#if it does not already exist, create a directory where the optimized dose will be saved
        os.mkdir(data_supradir+current+'/MRI/baseline/Prescribed_dose')
    out_dir = data_supradir+current+'/MRI/baseline/Prescribed_dose'
     
    print('Optimizing the dose prescription for '+current) 
              
    cell_map = sitk.ReadImage(data_supradir+current+'/MRI/baseline/micro/cellApp_CORRECTED_CTV.nii') #Read cellular density map. Units: um-3
    cell_map = cell_map*(10**8) #Convert cell map to Units: 10^4 cm-3
    
    cell_dir = data_supradir+current+'/MRI/baseline/micro/' #Path to dir with original data
    GTV_full = sitk.ReadImage(cell_dir+'GTV_inDWI_ITK.nii') #Read GTV ROI
                       
    GTV = sitk.ReadImage(cell_dir+'gtv_voxel_ok_3D.nii')
    GTV = sitk.Cast(GTV, sitk.sitkUInt8)
    GTV_outliers = GTV_full-GTV #Generate a mask of the outliers voxels by subtracting the voxels that have valid values of ADC from the CTV mask
    sitk.WriteImage(GTV_outliers, cell_dir+'GTV_outliers.nii')
    
    cell_map = generate_mask(cell_map, GTV) #Mask the cellularity map only in voxels of CTV that are valid for ADC
    cell_map_nda = sitk.GetArrayFromImage(cell_map) #Convert cell map image into a numpy array  
    
    #Find indexes of voxels that have values of cellularity
    result = np.where(cell_map_nda>0)
    listOfCoordinates= list(zip(result[0], result[1], result[2])) #Save indexes in a list of coordinates 

    orig_shape = sitk.GetArrayFromImage(cell_map).shape #Save the shape of the original cell map image
    voxel_volume = (cell_map.GetSpacing()[0]*cell_map.GetSpacing()[1]*cell_map.GetSpacing()[2])*10**(-3) # volume of each voxel in cm-3
    
    #Define input parameters for TCP modelling: number of cells per voxel (N0) and alpha and beta values per voxel    

    cell_dens = allVoxInt(cell_map, GTV) #Define array of values with cellular density>0
    N0 = cell_dens*voxel_volume #Multiply the values of cellular density per voxel by the volume of one voxel to obtain the number of cells per voxel

#%%Scale values of dose inside the GTV linearly with Cellularity according to: Di = Dmin + [(Dmax-Dmin)*(Celli-Cellmin)]/(Cellmax-Cellmin)
    
    Dmin = 74 #Gy minimum dose delivered to voxels in GTV with lowest cellularity
    Dmax = 81 #Gy maximum dose delivered to voxels in GTV with highest cellularity
    
    dose = linear_Dose_N0_mapping(N0, Dmin, Dmax)

#%% Saving the optimised dose values into a prescription dose image
    
    dose_nda = np.ndarray(shape=orig_shape) #Generate a 3D array of the same shape of the cell density map
    
    for d,i in zip(list(dose),listOfCoordinates): #Assign the calculated optimised dose to each voxel, based on the index of the corresponding cellularity voxel
        dose_nda[i]=d

    dose_img = sitk.GetImageFromArray(dose_nda) #Convert the numpy array into an image 
    dose_img.CopyInformation(cell_map) #Copy the headers of the cellularity image onto the headers of the optimised dose image
    
    # sitk.WriteImage(dose_img, out_dir+'/Dose_painted_GTV.nii') #Save the prescribed optimised dose image as a nifti file 
    
#%%Set values within GTV with cellularity untrusted to receive Dmin = 74 Gy
    
    dose_img_final = set_mask_value(dose_img, GTV_outliers, Dmin)
    sitk.WriteImage(dose_img_final, out_dir+'/Dose_painted_GTV_final.nii')
    
