## SBC Proton Dose Painting - Patients data pre-processing

### Author: Caterina Brighi

Pre-processing of the SBC patients data.
For a correct preprocessing pipeline, run scripts following the tasks order Task 1 - Task 7. 

#### Task 1: Reformatting RTSTRUCT dicom data for compatibility with RayStation: Plastimatch convert RTSTRUCT.txt

*    [x]  Reformatting RTSTRUCT dicom data for ensuring compatibility with RayStation when loading RTSTRUCT file with CT file. Run this in the terminal command line.
         
#### Task 2: Data conversion: dicom2nifti.py

*    [x]  Convert MRI, RTDOSE, RTSTRUCT and CT DICOM data to NIFTI format

#### Task 3: Coregister RT data to DWI MRI space: Transform_Targets_CTinDWIspace.py

*    [x]  Register CT image, GTV and CTV to T2 MRI image space, and then to DWI image space, by means applying given transformation matrices

#### Task 4: Patient's specific dose painting prescriptions derivation: Dpresc_fromCell_SBCproton.py

*    [x]  Generate a prescription dose map in DWI space, by linearly scaling the dose between 74-81 Gy with per-voxel cellularity in the GTV

#### Task 5: Patient's specific baseline boost dose prescriptions derivation: Generate_Boost_plans_prescriptions_SBCproton.py

*    [x]  Generate a prescription dose map in CT space, by assigning 74 Gy within the entire CTV

#### Task 6: Patient's specific inverse dose painting prescriptions derivation: Inverse_DosePrescriptions_SBCproton.py

*    [x]  Register the prescription dose map from DWI to CT space by applying the inverse transforms from Task 2
*    [x]  Generate inverse dose prescription maps in the GTV as Dinv = Dmax - Dpresc, where Dmax = maximum dose estimated in GTV for each patient
*    [x]  Derive CTV high dose (HD) margins and assign 74 Gy to all voxels in HD margins
*    [x]  Obtain a Dose prescription and an inverse dose prescription map in nifti format
*    [x]  Convert the inverse dose prescription nifti file into a RTdose dicom file

#### Task 7: Patient's specific inverse baseline boost dose prescriptions derivation: Inverse_DosePrescriptions_BoostPlans_SBCproton.py

*    [x]  Generate inverse dose prescription maps in the CTV as Dinv = Dmax - Dpresc, where Dmax = maximum dose estimated in CTV for each patient (i.e. 81 Gy) and Dpresc = 74 Gy
*    [x]  Derive CTV high dose (HD) margins and assign 74 Gy to all voxels in HD margins
*    [x]  Obtain an inverse dose prescription map in nifti format
*    [x]  Convert the inverse dose prescription nifti file into a RTdose dicom file

