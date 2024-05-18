
function ROIinfo = createROIinfoStruct(date)

% Creates structure relating ROI number, ADC, T2 for imaging on given date

switch date

    case '20240226'

        % Define ROI folder
        ROIfolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\")  num2str(date)  '\ROIs'];

        ROInums = {1,2, 3, 4, 8};
        ROInames = {'ROI1', 'ROI2', 'ROI3', 'ROI4', 'ROI8'};
        PVPconcs = {50, 40, 30, 20, 1};
        ADCs = {0.00012, 0.000227, 0.000382, 0.000581, 0.00101};
        T2s = {126.8, 235.3, 377.1, 539.1, 825.1};


    case '20240304'

        ROIfolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\")  num2str(date)  '\ROIs'];

        ROInums = {1,2, 3, 4, 5, 6, 7, 8};
        ROInames = {'ROI1', 'ROI2', 'ROI3', 'ROI4', 'ROI5', 'ROI6', 'ROI7', 'ROI8'};
        PVPconcs = {50, 40, 30, 20, 10, 5, 2.5, 1};
        ADCs = {0.00012, 0.000227, 0.000382, 0.000581, 0.000833, 0.000977, 0.001032, 0.00101};
        T2s = {126.8, 235.3, 377.1, 539.1, 693.6, 771.9, 807.3, 825.1};


    case '20240319'

        ROIfolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\")  num2str(date)  '\ROIs'];

        ROInums = {1, 2, 3, 4, 5, 6, 7, 8};
        ROInames = {'ROI1', 'ROI2', 'ROI3', 'ROI4', 'ROI5', 'ROI6', 'ROI7', 'ROI8'};
        PVPconcs = {50, 40, 30, 20, 10, 5, 2.5, 1};
        ADCs = {0.00012, 0.000227, 0.000382, 0.000581, 0.000833, 0.000977, 0.001032, 0.00101};
        T2s = {126.8, 235.3, 377.1, 539.1, 693.6, 771.9, 807.3, 825.1};
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
    ROIinfo(ROIindx).T2 = T2s{ROIindx};
    


end











end