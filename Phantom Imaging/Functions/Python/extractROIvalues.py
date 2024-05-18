# Python function to extract ROIs

import numpy as np
import SimpleITK as sitk
import sys
import glob
import os
from scipy.io import savemat

# Import DICOM from imgtools
sys.path.insert(0, r"C:\Users\adam\OneDrive - University College London\UCL PhD\Image-Processing\DICOM\Python")
import DICOM # type: ignore


def extractROIvalues(
    imagefname,
    ROIname,
    ROIftype = 'img',
    ROIfolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\ROIs",
    imagefolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\Images",
    savetype = 'mat',
    savefolder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\ROI values",
):
    
    '''
    Function to read ROI values from saved (and normalised) DWI image
    
    '''
    
    # Create dcm object for image
    DWIdcm = DICOM.MRIDICOM(imagefname)
    
    # Read series description
    SeriesDescription = DWIdcm.DICOM_df['Series Description'].values[0]
    
    # Set folder name with normalised diffusion images
    DWIimagefolder = f'{imagefolder}/{SeriesDescription}'
    
    # Find files and image names in folder
    imagefnames = glob.glob(f'{DWIimagefolder}/*.npy')
    imagenames = [os.path.split(imagefname)[1][:-4] for imagefname in imagefnames]
    print(DWIimagefolder)
    print(imagenames)
    # Load ROI
    print('yesss')
    ROIarray = sitk.GetArrayFromImage(sitk.ReadImage(f'{ROIfolder}/{SeriesDescription}/{ROIname}.{ROIftype}'))
    ROIbool = (ROIarray == 1)
    
    # For each image, extract values in ROI and save
    for imagename, imagefname in zip(imagenames, imagefnames):
        
        # Load normalised image 
        normimage = np.load(imagefname)
        
        # Extract ROI values
        ROIvals = normimage[ROIbool]
        
        # Save ROI values
        if savetype == 'mat':
            
            # Make directory
            folder = f'{savefolder}/{SeriesDescription}/mat/{ROIname}'
            try:
                os.makedirs(folder)
            except:
                None
            print('ok')
            # save 
            savemat(
                f'{folder}/{imagename}.mat',
                dict(zip(
                    ['ROIvals'], [ROIvals]
                ))
            )
            
            
                
                
                
                
                
                
                
                
                
                
                
    
extractROIvalues(
    imagefname, # type: ignore
    ROIname, # type: ignore
    ROIftype, # type: ignore
    ROIfolder, # type: ignore
    imagefolder, # type: ignore
    savetype, # type: ignore
    savefolder, # type: ignore
)