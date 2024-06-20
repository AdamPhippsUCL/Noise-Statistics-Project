import numpy as np
import SimpleITK as sitk
import sys
import glob
import os
from scipy.io import savemat

# Import DICOM from imgtools
sys.path.insert(0, r"C:\Users\adam\OneDrive - University College London\UCL PhD\Image-Processing\DICOM\Python")
import DICOM # type: ignore

# Define ADC filename
ADCfname = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Imaging Data\20230601_QIBAPHANTOM\DICOM\IM_0017"

# Create dcm object
ADCdcm = DICOM.MRIDICOM(ADCfname)

# Load ADC image
ADC = (ADCdcm.constructImageArray(VoxelType = 'DV'))*(1e-3)

# Load series description
name = ADCdcm.DICOM_df['Series Description'].values[0]
print(name)

# Save as mha
folder = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs\20230601_QIBAPHANTOM\ADC"
sitk.WriteImage(sitk.GetImageFromArray(ADC), f'{folder}/{name}.mha')

