% MATLAB script to save normalised DWI images

date = '20240319';

% Define DICOM folder
DICOMfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Imaging Data\20240319_AdamPhantom3\DICOM\DICOM";
        
% Define cell array of filenames
Imagecodes = {...
   'IM_0072',...
   'IM_0078'};

% Define multiframe bool
multiframe = true;

% Define save data type
datatype = 'mat';
% Define output folder 
outputfolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\") date char("\Images")];

% Create cell array of image fnames
Imagefnames = {};
for indx = 1:length(Imagecodes)
    Imagefnames{indx} = [char(DICOMfolder) '/' char(Imagecodes{indx})];
end

% 
% %% HARDCODING IMAGE FNAMES
% Imagefnames = {
%     "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Imaging Data\20240307_Patient1\DICOM\MSHP1 [MSHP1]\20000101 000000\Series 801 [MR - b1 NSA1DTI6]\1.3.6.1.4.1.5962.99.1.430794585.115837193.1709827811161.277.0.dcm",...
%     "C:\Users\adam\OneDrive - University College London\UCLPHD~1\PHDYEA~1\Projects\VERDIC~2\IMAGIN~1\202403~2\DICOM\MSHP1_~1\200001~1\SE5428~1\136141~2.DCM"...
% };
% 


%% Call Python script for each image fname
pyscript = "C:\Users\adam\OneDrive - University College London\UCL PhD\Image-Processing\DICOM\Python\normaliseDiffusionImage.py";

for indx = 1:length(Imagefnames)
    imagefname = Imagefnames{indx};
    pyrunfile( ...
        pyscript, ...
        imagefname = imagefname, ...
        multiframe = multiframe, ...
        returnimages = false,...
        saveimages = true,...
        datatype = datatype,...
        outputfolder = outputfolder...
        );
end