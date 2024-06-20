
function ROIinfo = createROIinfoStruct(date, ROIfolder)

% Creates structure relating ROI number, ADC, T2 for imaging on given date

switch date

    % case '20240226'
    % 
    %     ROInums = {1,2, 3, 4, 8};
    %     ROInames = {'ROI1', 'ROI2', 'ROI3', 'ROI4', 'ROI8'};
    %     PVPconcs = {50, 40, 30, 20, 1};
    %     ADCs = {0.00012, 0.000227, 0.000382, 0.000581, 0.00101};
    %     T2s = {126.8, 235.3, 377.1, 539.1, 825.1};
    % 
    % 
    % case '20240304'
    % 
    %     ROInums = {1,2, 3, 4, 5, 6, 7, 8};
    %     ROInames = {'ROI1', 'ROI2', 'ROI3', 'ROI4', 'ROI5', 'ROI6', 'ROI7', 'ROI8'};
    %     PVPconcs = {50, 40, 30, 20, 10, 5, 2.5, 1};
    %     ADCs = {0.00012, 0.000227, 0.000382, 0.000581, 0.000833, 0.000977, 0.001032, 0.00110};
    %     T2s = {126.8, 235.3, 377.1, 539.1, 693.6, 771.9, 807.3, 825.1};


    case '20240319'

        ROInums = {1, 2, 3, 4, 5, 6, 7, 8};
        ROInames = {'ROI1', 'ROI2', 'ROI3', 'ROI4', 'ROI5', 'ROI6', 'ROI7', 'ROI8'};
        PVPconcs = {50, 40, 30, 20, 10, 5, 2.5, 1};
        ADCs = {0.12e-3, 0.25e-3, 0.382e-3, 0.581e-3, 0.833e-3, 0.977e-3, 1.032e-3, 1.101e-3};
        % ADCerrs = {0.006e-3, 0.007e-3, 0.01e-3, 0.015e-3, 0.023e-3, 0.023e-3, 0.025e-3, 0.027e-3};
        T2s = {126.8, 235.3, 377.1, 539.1, 693.6, 771.9, 807.3, 825.1};
        % T2errs = {2.2, 4.1, 6.6, 9.4, 12.1, 13.5, 14.1, 14.4 };

        % % MY MEASUREMENTS
        % ADCs = {1.196e-4, 2.348e-4, 3.76e-4, 5.805e-4, 8.06e-4, 9.47e-4, 1.01e-3, 1.095e-3};
        % ADClbs = {1e-4, 2.229e-4, 3.554e-4, 5.618e-4, 7.77e-4, 9.18e-4, 9.73e-4, 1.075e-3};
        % ADCubs = {1.4e-4, 2.47e-4, 3.94e-4, 6.038e-4,  8.34e-4, 9.88e-4, 1.075e-3, 1.12e-3};
        % 
        % T2s = {132.5, 254, 426, 590, 716, 762, 806, 881};
        % T2lbs = {123.9, 243, 413, 560, 681, 728, 760, 802 };
        % T2ubs = {141.0, 266, 440, 620, 750, 801, 860, 960};

end


%% Create struct

ROIinfo = struct();


NROI = length(ROInums);

for ROIindx = 1:NROI
    
    ROIinfo(ROIindx).ROIfolder = ROIfolder;
    ROIinfo(ROIindx).ROInum = ROInums{ROIindx};
    ROIinfo(ROIindx).ROIname = ROInames{ROIindx};
    ROIinfo(ROIindx).PVPconc = PVPconcs{ROIindx};
    ROIinfo(ROIindx).ADC = ADCs{ROIindx};
    % ROIinfo(ROIindx).ADClb = ADClbs{ROIindx};
    % ROIinfo(ROIindx).ADCub = ADCubs{ROIindx};
    ROIinfo(ROIindx).T2 = T2s{ROIindx};
    % ROIinfo(ROIindx).T2lb = T2lbs{ROIindx};    
    % ROIinfo(ROIindx).T2ub = T2ubs{ROIindx};  

end











end