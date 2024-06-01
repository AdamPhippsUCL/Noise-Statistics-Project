function ImageSettings =  createImageSettingStruct(date, experiment)

arguments

    date % Date of phantom imaging session
    experiment % experiment number (details below)

end

%% DATA

% ====== DATE: 20240226

% == Experiment 1: b90 vary NSA

% Image numbers:

% IM_0003 : b90_vx1.3
% IM_0006 : b90 4x
% IM_0009 : b90 2x
% IM_0012 : b90 0p5x
% IM_0015 : b90 1p5x

% == Experiment 2: b500 vary NSA

% Image numbers

% IM_0018 : b500_NSA4
% IM_0021 : b500_NSA4x4
% IM_0024 : b500_NSA4x2
% IM_0027 : b500_NSA4x0p5


% == Experiment 3: VERDICT protocol (varying b and TE)

% Image numbers

% IM_0003 : b90_vx1.3
% IM_0018 : b500_NSA4
% IM_0030 : b1500_NSA4
% IM_0033 : b2000_NSA4
% IM_0036 : b3000_NSA4


% Experiment 4 and 5: DTI for sigma0 measurement


% ===== DATE: 20240304

% == Experiment 1

% IM_0003 : b2000_highb3
% IM_0006 : b2000_highb9
% IM_0009 : b2000_highb6
% IM_0012 : b2000_highb18





% ===== DATE: 20240319

% == Experiment 1: varying TE/b settings

% IM_0006 : b2000_Ex1
% IM_0009 : b2500_Ex1
% IM_0015 : b3000_Ex1
% IM_0018 : b4000_Ex1
% IM_0021 : b5000_Ex1
% IM_0024 : b2750_Ex1
% IM_0027 : b1500_Ex1
% IM_0030 : b3500_Ex1
% IM_0033 : b2100_Ex1
% IM_0036 : b2400_Ex1

% == Experiment 2: vary NSA

% IM_0009 : b2500_Ex1
% IM_0039 : b2500_NSA2_Ex2
% IM_0042 : b2500_NSA4_Ex2
% IM_0045 : b2500_NSA6_Ex2
% IM_0048 : b2500_NSA8_Ex2

% == Experiment 3: vary Avg. high b

% IM_0009 : b2500_Ex1
% IM_0051 : b2500_Avg1_Ex3
% IM_0054 : b2500_Avg2_Ex3
% IM_0057 : b2500_Avg4_Ex3
% IM_0060 : b2500_Avg6_Ex3
% IM_0063 : b2500_Avg8_Ex3



%% Define scan data

switch date

    case '20240226'

        dicomfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Imaging Data\20240226_noise\DICOMSall\DICOM";

        switch experiment

            case 1

                image_nums = {'IM_0003', 'IM_0006', 'IM_0009', 'IM_0012', 'IM_0015'};
                bvals = {90, 90, 90, 90, 90};
                TEs = {55.2, 55.2, 55.2, 55.2, 55.2};
                NSAs = {4, 16, 8, 2, 6};
                Nav_ratios = {1,1,1,1,1};

            case 2

                image_nums = {'IM_0018', 'IM_0021', 'IM_0024', 'IM_0027'};
                bvals = {500, 500, 500, 500};
                TEs = {68, 68, 68, 68};
                NSAs = {4, 16, 8, 2};
                Nav_ratios = {2,2,2,2};            

            case 3

                image_nums = {'IM_0003', 'IM_0018', 'IM_0030', 'IM_0033', 'IM_0036'};
                bvals = {90, 500, 1500, 2000, 3000};
                TEs = {55.2, 68, 94, 76, 87};
                NSAs = {4, 4, 4, 4, 4};
                Nav_ratios = {1, 2, 3, 3, 3};                      


        end


    case '20240304'

        dicomfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Imaging Data\20240304_Vnoise\DICOM";

        switch experiment

            case 1

                image_nums = {'IM_0003', 'IM_0009', 'IM_0006' ,'IM_0012'};
                bvals = {2000,2000,2000,2000};
                TEs = {76, 77, 76, 77};
                NSAs = {4,4,4,4};
                Nav_ratios = {3, 6, 9, 18};
    

            case 2 

                image_nums = {'IM_0015', 'IM_0018', 'IM_0021', 'IM_0024', 'IM_0027', 'IM_0030', 'IM_0033', 'IM_0036', 'IM_0039', 'IM_0042', 'IM_0045'};
                bvals =      {90,         500,       1500,      2000,      3000,      800,       1200,      1600,      1800,      2400,      1600};
                TEs =        {55,         68,        94,        76,        87,        80,        70,        80,        83,        91,        80};
                NSAs =       {4,          4,         4,         4,         4,         4,         4,         4,         4,         4,         4};
                Nav_ratios = {3,          6,         9,         9,         9,         6,         9,         9,         9,         9,         3};
                

        end


    case '20240319'
        
        dicomfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Imaging Data\20240319_AdamPhantom3\DICOM\DICOM";
        
        switch experiment

            case 1
                
                
                image_nums = {'IM_0006', 'IM_0009', 'IM_0015',  'IM_0024', 'IM_0027', 'IM_0030', 'IM_0033', 'IM_0036', 'IM_0018', 'IM_0021'};
                bvals = {2000, 2500, 3000, 2750, 1500, 3500, 2100, 2400, 4000, 5000};
                TEs = {200, 300, 400, 450, 450, 150, 100, 234, 250, 325};
                NSAs = {1,1,1,1,1,1,1,1,1,1};
                Nav_ratios = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
            

            case 2

                image_nums = {'IM_0051', 'IM_0054', 'IM_0057', 'IM_0060', 'IM_0063'};
                bvals = { 2500, 2500, 2500, 2500, 2500};
                TEs = { 300, 300, 300, 300, 300};
                NSAs = {2,2,2,2,2};
                Nav_ratios = { 1*3, 2*3, 4*3, 6*3, 8*3};

            case 3

                image_nums = {'IM_0009', 'IM_0039', 'IM_0042', 'IM_0045', 'IM_0048'};
                bvals = {2500, 2500, 2500, 2500, 2500};
                TEs = {300, 300, 300, 300, 300};
                NSAs = {1, 2, 4, 6, 8};
                Nav_ratios =  {3, 3, 3, 3, 3};

        end

end


%% Create structure

ImageSettings = struct();

Nimage = length(image_nums);

for ImIndx = 1:Nimage

    ImageSettings(ImIndx).imagefname = [char(dicomfolder) '/' image_nums{ImIndx}];
    ImageSettings(ImIndx).bval = bvals{ImIndx};
    ImageSettings(ImIndx).TE = TEs{ImIndx};
    ImageSettings(ImIndx).NSA = NSAs{ImIndx};
    ImageSettings(ImIndx).Nav_ratio = Nav_ratios{ImIndx};

end








end