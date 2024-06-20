import numpy as np
import SimpleITK as sitk
import sys
import glob
import os
from scipy.io import savemat

# Import DICOM from imgtools
sys.path.insert(0, r"C:\Users\adam\OneDrive - University College London\UCL PhD\Image-Processing\DICOM\Python")
import DICOM # type: ignore


# Define LWI fname
LWIfname = r"C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Imaging Data\20240304_Vnoise\DICOM\IM_0079"

# Create DCM object
LWIdcm = DICOM.MRIDICOM(LWIfname)

# Image array
LWI = LWIdcm.constructImageArray()

# Save as mha
sitk.WriteImage(sitk.GetImageFromArray(LWI), rf'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs\20240304_Vnoise\Images\LWI.mha')

# TEs
TEs = LWIdcm.DICOM_df['Effective Echo Time'].values
TEvals = list(set(TEs))
print(TEvals)
# # Separate images with different TEs
# for TEval in TEvals:

#     TEimg = LWI[TEs == TEval]

#     # Save as mha
#     sitk.WriteImage(sitk.GetImageFromArray(TEimg), rf'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs\20240304_Vnoise\Images\TE_{TEval}.mha')




# T2 calculation
T2results = np.zeros((8, 3))

for ROIindx in range(1,9):

    # Load ROI
    ROI = sitk.GetArrayFromImage(sitk.ReadImage(rf'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs\20240304_Vnoise\ROIs\ROI{ROIindx}.img'))
    ROI = (ROI==1)

    # Load values from different echo times
    for indx, TEval in enumerate(TEvals):

        thisTEimg = LWI[TEs == TEval]

        thisTE_ROIvals = (thisTEimg[ROI]).flatten()

        if indx == 0:
            ROIvals = np.zeros( (8,len(thisTE_ROIvals)))

        ROIvals[indx,:] = thisTE_ROIvals


    # Take log of data
    logROIvals = (np.log(ROIvals))

    # Straight line fit
    coeffs, res = np.polyfit(x=np.transpose(TEvals), y=logROIvals, deg = 1)

    T2s = -1/coeffs

    T2results[ROIindx-1,0] = np.mean(T2s)
    T2results[ROIindx-1,1] = np.percentile(T2s,5)
    T2results[ROIindx-1,2] = np.percentile(T2s, 95)

print(T2results)