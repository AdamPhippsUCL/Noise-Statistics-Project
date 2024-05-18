% MATLAB script to apply fitting to histograms

% Imaging date
date = 20240319;

% Image fname
imagefname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Imaging Data\20240319_AdamPhantom3\DICOM\DICOM\IM_0039";
% Define scan info
bval = 2500;
TE = 300;
Nav_ratio = 3;

% Distribution type
disttype = 'Rice';


% Define ROI folder
ROIfolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\")  num2str(date)  '\ROIs'];

% Define ROI number 
ROIname = 'ROI1';

% ROI info
ADC = 0.12e-3;
T2 = 126;


%% Check if all data exists

% Define series description
dinfo = dicominfo(imagefname);
SeriesDescription = dinfo.SeriesDescription;

% Define image folder
imagefolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\")  num2str(date)  '\Images\numpy'];

% Check if normalised images are saved, if not: go and save image!
if ~exist([imagefolder '/' SeriesDescription ], 'dir')
    error('Normalised image not saved!')
end

% Check if ROI exists, if not: define one!
if ~exist([ROIfolder '/' SeriesDescription '/' ROIname '.img'], 'file')
    error(['ROI: ' ROIname ' not saved!'])
end


%% Extract ROI values
ROIvalsfolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\")  num2str(date)  '\ROI values'];

% % Check if ROI values saved
% if ~exist([ROIvalsfolder '/' SeriesDescription '/mat/' ROIname ], 'dir')
%     extractROIvalues(imagefname, ROIname);
% end
extractROIvalues( ...
    imagefname, ...
    ROIname,...
    imagefolder = imagefolder,...
    ROIfolder = ROIfolder,...
    savefolder = ROIvalsfolder);

% Define diffusion orientation
DiffDirec = 0;

% Load ROI values
load([ROIvalsfolder '/' SeriesDescription '/mat/' ROIname '/bval' num2str(bval) '.0_DiffDirec' num2str(DiffDirec) '.mat'])


%% Make histogram

% Define bin edges and centers
binmin = 0;
binmax = 2;
nbin = 150;
binedges = linspace(binmin, binmax, nbin+1);
bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
binspacing = bincenters(2)-bincenters(1);

figure;
H = histogram(ROIvals, binedges);
hold on;
counts = H.Values;


%% Fit to histogram

bval = bval;
% Predicted diffusion signal fraction
fd =exp(-ADC*bval);

% Initial guess
beta0guess = [std(ROIvals)/sqrt(2), T2, exp(-ADC*bval)];

% lowerbound
lb = [0.005, 0.8*T2, 0.8*exp(-ADC*bval)];
ub = [0.1, 1.2*T2, 1.2*exp(-ADC*bval)];


% Set parameters to NaN if not fixed in fitting
T2 = NaN;
fd = NaN;



[coeffs, resnorm] = fitDistToHist( ...
    counts, ...
    bincenters, ...
    T2 = T2, ...
    fd = fd, ...
    TE = TE, ...
    Nav_ratio=Nav_ratio,...
    disttype = disttype, ...
    beta0guess = beta0guess,...
    lb = lb,...
    ub=ub);


sigma0 = coeffs(1);
T2_predict = coeffs(2);
fd_predict = coeffs(3);

%% Generate pdf

b0signal = exp(-TE/T2_predict);
bsignal = fd_predict*b0signal;

switch disttype
    case 'Ratio'
        [dist, signals] = RatioDistRician(b0signal, bsignal, sigma0, Nav_ratio=Nav_ratio);
    case 'Rice'
        [dist, signals] = RiceDist(b0signal, bsignal, sigma0, Nav_ratio=Nav_ratio);
end

hold on
plot(signals, dist*binspacing*sum(counts), '-*');


predicted_sigma0 = sigma0
predicted_T2 = T2_predict
predicted_ADC = (-1/bval)*log(fd_predict)
resnorm

