# -*- coding: utf-8 -*-
"""
Created on Thu May 13 11:58:57 2021

@author: cbri3325
"""



import math
import re
import struct
import numpy as np
import os
import time
import pydicom
import SimpleITK as sitk



GTransferSyntaxUID = "1.2.840.10008.1.2"
GImplementationClassUID = "1.2.826.0.1.3680043.8.498.75006884747854523615841001"

RTDOSEModality = "RTDOSE"
RTPLANModality = "RTPLAN"
RTSTRUCTModality = "RTSTRUCT"

RTDoseSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.2"
RTStructSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3"
RTPlanSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.5"

Manufacturer = "RayStation"

# dose = sitk.ReadImage('C:/Users/cbri3325/Dropbox (Sydney Uni)/Caterina Brighi/Data/test/dose.nii')
# reference_dcm = 'C:/Users/cbri3325/Dropbox (Sydney Uni)/Caterina Brighi/Data/test/RD.zzSPARK_PAT05.Dose_PLAN.dcm'
# Inv_DP_origCT_inPTV
# dcm_ref

def convert_nifti_to_dicom_RTdosefile(image, reference_dcm, output_directory=".", out_filename="dose.dcm"):
    """Converts a Nifti image to a Dicom RT dose file
    Args:
        image (sitk.Image): A SimpleITK image object to convert
        reference_dcm (str): A directory path containing a reference Dicom RT dose file to use
        output_directory (str, optional): The directory in which to place the generated Dicom
                                          files. Defaults to ".".
    """

    # Make the output directory if it doesn't already exist
    os.makedirs(output_directory, exist_ok=True)

    # Read in the reference Dicom RT dose file
    dicom_object = pydicom.read_file(reference_dcm, force=True)
    # raw_dose_image = sitk.ReadImage(reference_dcm, sitk.sitkFloat32)
    # dose_grid_scaling = dicom_object.DoseGridScaling
    # scaled_dose_image = raw_dose_image * dose_grid_scaling

    
    # Generate some new UIDs
    doseInstanceUID = pydicom.uid.generate_uid()
    
    # Populate required values for file meta information
    file_meta = pydicom.dataset.Dataset()
    file_meta.MediaStorageSOPClassUID = RTDoseSOPClassUID
    file_meta.TransferSyntaxUID = GTransferSyntaxUID
    file_meta.MediaStorageSOPInstanceUID = doseInstanceUID
    file_meta.ImplementationClassUID = GImplementationClassUID
    
    # Create the pydicom.dataset.FileDataset instance (initially no data elements, but file_meta supplied)
    RDfilename = f"RD.{file_meta.MediaStorageSOPInstanceUID}.dcm"
    ds = pydicom.dataset.FileDataset(
        RDfilename, {}, file_meta=file_meta, preamble=b"\x00" * 128
    )
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")

    ds.SOPClassUID = RTDoseSOPClassUID  # RT Dose Storage
    ds.SOPInstanceUID = doseInstanceUID
    ds.StudyDate = dicom_object.StudyDate
    ds.StudyTime = dicom_object.StudyTime
    ds.AccessionNumber = ""
    ds.Modality = "RTDOSE"
    ds.Manufacturer = dicom_object.Manufacturer
    ds.ReferringPhysicianName = dicom_object.ReferringPhysicianName
    # ds.StationName = dicom_object.StationName
    # ds.StudyDescription = dicom_object.StudyDescription
    ds.SeriesDescription = 'Inverse dose prescription'
    ds.ManufacturerModelName = dicom_object.ManufacturerModelName
    ds.PatientName = dicom_object.PatientName
    ds.PatientID = dicom_object.PatientID
    ds.PatientBirthDate = dicom_object.PatientBirthDate
    ds.PatientSex = dicom_object.PatientSex
    ds.SliceThickness = int(image.GetSpacing()[2])
    # ds.SliceThickness = dicom_object.SliceThickness
    # ds.DeviceSerialNumber = dicom_object.DeviceSerialNumber
    # ds.SoftwareVersion = dicom_object.SoftwareVersion
    ds.StudyInstanceUID = dicom_object.StudyInstanceUID
    ds.SeriesInstanceUID = dicom_object.SeriesInstanceUID
    ds.StudyID = dicom_object.StudyID
    ds.SeriesNumber = dicom_object.SeriesNumber
    ds.InstanceNumber = dicom_object.InstanceNumber
    ds.ImagePositionPatient = list(image.GetOrigin()) 
    ds.ImageOrientationPatient = list(image.GetDirection()[0:6])
    # ds.ImagePositionPatient = dicom_object.ImagePositionPatient
    # ds.ImageOrientationPatient = dicom_object.ImageOrientationPatient
    ds.FrameOfReferenceUID = dicom_object.FrameOfReferenceUID
    ds.PositionReferenceIndicator = dicom_object.PositionReferenceIndicator
    ds.SamplesPerPixel = image.GetNumberOfComponentsPerPixel()
    # ds.SamplesPerPixel = dicom_object.SamplesPerPixel
    ds.PhotometricInterpretation = dicom_object.PhotometricInterpretation
    # ds.NumberOfFrames = dicom_object.NumberOfFrames
    ds.NumberOfFrames = int(image.GetSize()[2])
    ds.FrameIncrementPointer = dicom_object.FrameIncrementPointer
    ds.Rows = int(image.GetSize()[1])
    ds.Columns = int(image.GetSize()[0])
    ds.PixelSpacing = list(image.GetSpacing()[0:2])
    # ds.Rows = dicom_object.Rows
    # ds.Columns = dicom_object.Columns
    # ds.PixelSpacing = dicom_object.PixelSpacing
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = dicom_object.PixelRepresentation
    ds.DoseUnits = dicom_object.DoseUnits #'GY'
    ds.DoseType = dicom_object.DoseType #'PHYSICAL'
    ds.DoseSummationType = dicom_object.DoseSummationType #'PLAN'
    slice_thickness = image.GetSpacing()[2]
    ds.GridFrameOffsetVector = [ds.ImagePositionPatient[2] + x*slice_thickness for x in range(ds.NumberOfFrames)]

    np_dose = sitk.GetArrayFromImage(image)
    max_dose_val = np_dose.max()

    ds.DoseGridScaling = max_dose_val/65536 #Divide maximum dose value of the image by 2^16
    print('Dose grid scaling factor: '+str(ds.DoseGridScaling))
    np_dose_scaled = np_dose/ds.DoseGridScaling
    np_dose_scaled = np_dose_scaled.astype(np.uint16)
    ds.TissueHeterogeneityCorrection = "IMAGE"
    # ds.ReferencedRTPlanSequence = dicom_object.ReferencedRTPlanSequence #need to link to RayStation frame of reference
    # ds.ReferencedSOPClassUID = dicom_object.ReferencedSOPClassUID
    
    #Need to copy intensity voxel values to pixels into dose image  
    ds.PixelData = np_dose_scaled.tobytes()

    # Save the RTDose Dicom File
    output_file = os.path.join(output_directory, out_filename)
    ds.save_as(output_file)

if __name__ == "__main__":
    convert_nifti_to_dicom_RTdosefile("./Data")