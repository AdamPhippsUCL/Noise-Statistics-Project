function extractROIvalues(imagefname, ROIname, opts)

arguments
    imagefname % filnemae of image DICOM
    ROIname % name of ROI

    % Options
    opts.ROIftype = 'img'
    opts.ROIfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\20240226\ROIs"
    opts.savetype = 'mat' % data type to save ROI values
    opts.savefolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\20240226\ROI values"
    opts.imagefolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\20240226\Images"
    opts.pyfile = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Phantom Imaging\Functions\Python\extractROIvalues.py"

end

% Function saves ROI values for each normalised diffusion image associated
% with image filename

pyrunfile( ...
    opts.pyfile, ...
    imagefname = imagefname,...
    ROIname = ROIname,...
    ROIfolder = opts.ROIfolder,...
    ROIftype = opts.ROIftype,...
    savetype = opts.savetype,...
    savefolder = opts.savefolder,...
    imagefolder = opts.imagefolder)




end