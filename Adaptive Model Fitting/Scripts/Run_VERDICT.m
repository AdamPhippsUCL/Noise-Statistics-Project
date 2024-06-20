% MATLAB script to run VERDICT processing on specified patients

%% Define study path (path to folder containing patients)

% Patient volunteers (Short VERDICT)
STUDY_path = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Imaging Data\Patient Volunteers\PE P";

%% Define patient numbers

% FOR VOLUNTEERS
PatNums = {"20240610"};

%% SAVED DATA

% Decide whether to used presaved MATLAB data to speed up processing
UseSavedData = true;

%% DEFINE FOLDERS

% Define output folder
OutputFolder= "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs";

% Define schemes folder
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Adaptive Model Fitting\Schemes";

% Define models folder
modelsfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Adaptive Model Fitting\MLP Models";

% Define python folder
pythonfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Adaptive Model Fitting\Functions\Python";

%% DEFINE VERDICT PROTOCOL

% === Model type
% modeltypes = { 'Original VERDICT', 'No VASC VERDICT'}; 
modeltypes = {'No VASC VERDICT'};

% === Scheme name
% schemenames = { 'Original Full', 'Original ex905003000'} ;
schemenames = {'Short Scheme v1'};

% === fitting technique
% fittingtechniques =   { 'AMICO', 'MLP'};
fittingtechniques = {'MLP'};

% === Noise used in MLP training
noisetype='Rice';
sigma0train = 0.025;
T2train = 10000;

% === ADC
calcADC = true;
% Max b value
vADCbmax = 1501;

%% Noise calibration

% Define fnames of dual echo images
echo1fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Imaging Data\Patient Volunteers\PE P\20240610\Series 1101 [MR - b1 NSA1LOW6te50P]\1.3.6.1.4.1.5962.99.1.740795050.145724320.1718727746218.570.0.dcm";
echo2fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Imaging Data\Patient Volunteers\PE P\20240610\Series 1301 [MR - b1 NSA1LOW6te125P]\1.3.6.1.4.1.5962.99.1.740795050.145724320.1718727746218.820.0.dcm";
echofnames = {echo1fname, echo2fname};

% RUN NOISE CALIBRATION HERE, T2 and sigma0 as inputs 


for protindx = 1:length(modeltypes)
    modeltype = modeltypes{protindx};
    schemename = schemenames{protindx};
    fittingtechnique = fittingtechniques{protindx};

    % Run VERDICT processing
    for patindx = 1:size(PatNums,1)
        PatNum = PatNums{patindx};
        disp(["----------->" PatNums(patindx)])
        VERDICT( ...
            PatNum, ...
            UseSavedData = UseSavedData,...
            modeltype = modeltype,...
            schemename = schemename,...
            fittingtechnique = fittingtechnique,...
            noisetype=noisetype,...
            sigma0train = sigma0train,...
            T2train=T2train,...
            STUDY_path=STUDY_path,...
            parent_folder=OutputFolder,...
            schemesfolder = schemesfolder,...
            modelsfolder = modelsfolder,...
            pythonfolder = pythonfolder,...
            calcADC = calcADC,...
            vADCbmax = vADCbmax...
            );
    end

end


