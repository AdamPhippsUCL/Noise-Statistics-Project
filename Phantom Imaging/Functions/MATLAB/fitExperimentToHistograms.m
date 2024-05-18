function FittingResults = fitExperimentToHistograms(date, experiment, opts)

arguments
    date % Date of phantom imaging session
    experiment % experiment number (details below)

    % OPTIONS
    opts.ROInums = [1,2,3,4,5,6,7,8] % Specify list of ROI numbers if not all analysed
    opts.disttype = 'Rice'
    opts.nbin= 125
    opts.paramtolerance = 0.1 % Fraction tolerance in fitted T2 and ADC values
    opts.sigma0max = 0.25
    opts.DiffDirec = 0 % Diffusion direcion
    opts.mainfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\"
end

%% Image and ROI details

ImageSettings = createImageSettingStruct(date, experiment);
Nimage = length(ImageSettings);

ROIinfo= createROIinfoStruct(date);

% SELECT ROIs desired!!!!
if ~strcmp(opts.ROInums, 'all')   
    for ROInum = [ROIinfo.ROInum]  
        if sum(ROInum == opts.ROInums) == 0
            nums = [ROIinfo.ROInum];
            bools = (nums == ROInum);
            ROIinfo(bools) = [];          
        end
    end
end

NROI = length(ROIinfo);


%% FIT HISTOGRAMS

FittingResults = struct();

% Loop over images and ROIs
for imageIndx = 1:Nimage

    % Get image settings
    imagefname = ImageSettings(imageIndx).imagefname;
    bval = ImageSettings(imageIndx).bval;
    TE = ImageSettings(imageIndx).TE;
    NSA = ImageSettings(imageIndx).NSA;
    Nav_ratio = ImageSettings(imageIndx).Nav_ratio;


    % Define series description
    dinfo = dicominfo(imagefname);
    SeriesDescription = dinfo.SeriesDescription;

    % Define main image folder
    mainimagefolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\")  num2str(date)  '\Images'];
        
    % Define numpy image folder
    imagefolder = [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics\")  num2str(date)  '\Images\numpy'];
    
    % Check if normalised images are saved, if not: go and save image!
    if ~exist([imagefolder '/' SeriesDescription ], 'dir')
        error(['Normalised image not saved: ' SeriesDescription])
    end   


    for ROIindx = 1:NROI

        % Get ROI settings
        ROIfolder = ROIinfo(ROIindx).ROIfolder;
        ROIname = ROIinfo(ROIindx).ROIname;
        ROInum = ROIinfo(ROIindx).ROInum;
        ADC = ROIinfo(ROIindx).ADC;
        T2 = ROIinfo(ROIindx).T2;

        % Check if ROI exists, if not: define one!
        if ~exist([ROIfolder '/' SeriesDescription '/' ROIname '.img'], 'file')
            error(['ROI: ' ROIname ' not saved for image: ' SeriesDescription])
        end



        %% Extract ROI values
        ROIvalsfolder = [char(opts.mainfolder) num2str(date)  '\ROI values'];
        
        % % Check if ROI values saved
        % if ~exist([ROIvalsfolder '/' SeriesDescription '/mat/' ROIname ], 'dir')
        %     extractROIvalues( ...
        %         imagefname, ...
        %         ROIname,...
        %         imagefolder = imagefolder,...
        %         ROIfolder = ROIfolder,...
        %         savefolder = ROIvalsfolder);
        % end

        extractROIvalues( ...
            imagefname, ...
            ROIname,...
            imagefolder = imagefolder,...
            ROIfolder = ROIfolder,...
            savefolder = ROIvalsfolder);

        % Define diffusion orientation
        DiffDirec = opts.DiffDirec;
        
        % Load ROI values
        load([ROIvalsfolder '/' SeriesDescription '/mat/' ROIname '/bval' num2str(bval) '.0_DiffDirec' num2str(DiffDirec) '.mat'])

        disp(std(ROIvals))

        %% Make histogram
        
        % Define bin edges and centers
        binmin = 0;
        binmax = 2;
        nbin = opts.nbin;
        binedges = linspace(binmin, binmax, nbin+1);
        bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
        binspacing = bincenters(2)-bincenters(1);
        
        f = figure;%('visible','off');
        H = histogram(ROIvals, binedges, FaceColor = [0 0.4470 0.741], FaceAlpha = 0.25);
        hold on;
        counts = H.Values;
        % close(f);


        % === Fit to histogram
      
        % Predicted diffusion signal fraction
        fd = exp(-ADC*bval);
        
        % Initial guess
        beta0guess = [std(ROIvals), T2, exp(-ADC*bval)];
        
        % lowerbound
        a = (1-opts.paramtolerance);
        b = (1+opts.paramtolerance);

        lb = [0.001, 0.9999*T2, exp(-b*ADC*bval)];
        if 0.999>exp(-a*ADC*bval)
            ub = [opts.sigma0max, 1.00001*T2,exp(-a*ADC*bval)];
        else
            ub = [opts.sigma0max, 1.00001*T2,0.9999];
        end
        
        % Set parameters to NaN if not fixed in fitting
        T2 = NaN;
        fd = NaN;       

        % Distribution coefficients
        [coeffs, resnorm] = fitDistToHist( ...
            counts, ...
            bincenters, ...
            T2 = T2, ...
            fd = fd, ...
            TE = TE, ...
            Nav_ratio=Nav_ratio,...
            disttype = opts.disttype, ...
            beta0guess = beta0guess,...
            lb = lb,...
            ub=ub);

        sigma0 = coeffs(1);
        T2fit= coeffs(2);
        fdfit = coeffs(3);

        ADCfit = (-1/bval)*log(fdfit);


        %% Plot fitted distribution
        b0signal = exp(-TE/T2fit);
        bsignal = fdfit*b0signal;
        dz = (binmax-binmin)/nbin;

        switch opts.disttype
            case 'Ratio'
                [dist, signals] = RatioDistRician(b0signal, bsignal, sigma0, Nav_ratio=Nav_ratio, zmin = binmin, zmax = binmax , dz = dz);
            case 'Rice'
                [dist, signals] = RiceDist(b0signal, bsignal, sigma0, Nav_ratio=Nav_ratio, zmin = binmin, zmax = binmax , dz = dz);
        end
        
        hold on
        
        xlabel('Normalised signal')
        ylabel('Counts')
        plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = opts.disttype, color = "#D95319");

        close(f);

        % Fill out fitting results structure
        FittingResults((imageIndx-1)*NROI + ROIindx).SeriesDescription = SeriesDescription;
        FittingResults((imageIndx-1)*NROI + ROIindx).ImageSettings = ImageSettings(imageIndx);
        FittingResults((imageIndx-1)*NROI + ROIindx).imageIndx = imageIndx;
        FittingResults((imageIndx-1)*NROI + ROIindx).ROIindx = ROIindx;
        FittingResults((imageIndx-1)*NROI + ROIindx).ROInum = ROInum;        
        FittingResults((imageIndx-1)*NROI + ROIindx).ROInum = ROIname;        
        FittingResults((imageIndx-1)*NROI + ROIindx).ROIinfo = ROIinfo(ROIindx);
        FittingResults((imageIndx-1)*NROI + ROIindx).sigma0 = sigma0;
        FittingResults((imageIndx-1)*NROI + ROIindx).T2fit = T2fit;
        FittingResults((imageIndx-1)*NROI + ROIindx).ADCfit = ADCfit;
        FittingResults((imageIndx-1)*NROI + ROIindx).resnorm = resnorm;
    end
end


%% Plot results

% Three figures: for sigma0, T2, and ADC fits respectively

% === Figure 1: sigma0 estimation

% Box plot for sigma0 estimation in each ROI

fig1 = figure;

ROInames  = [];
Sigma0s = zeros(NROI, Nimage);

for ROIindx = 1:NROI

    % Get ROI settings
    ROIname = ROIinfo(ROIindx).ROIname;
    ROInames = [ROInames string(ROIname)];
    ROInum = ROIinfo(ROIindx).ROInum;

    % ROI bools
    ROIbools = ([FittingResults.ROIindx] == ROIindx);

    % Extract sigma0 measurements
    Sigma0s(ROIindx, :) = [FittingResults(ROIbools).sigma0];
    scatter(ROIindx*ones(Nimage, 1), Sigma0s(ROIindx, :), '*')
    hold on

end

boxplot(transpose(Sigma0s))
xticks(linspace(1,NROI,NROI))
xticklabels(ROInames)
ylabel('Estimated \sigma_0 value')
ylim([0, opts.sigma0max])
grid on
title(opts.disttype)


% === Figure 2: T2 estimation

% Box plot for ADC estimation in each ROI

fig2 = figure;

ROInames  = [];
T2prcnterrs = zeros(NROI, Nimage);

for ROIindx = 1:NROI

    % Get ROI settings
    ROIname = ROIinfo(ROIindx).ROIname;
    ROInames = [ROInames string(ROIname)];
    ROInum = ROIinfo(ROIindx).ROInum;

    % ROI bools
    ROIbools = ([FittingResults.ROIindx] == ROIindx);

    % Extract T2 measurements error presentages
    T2prcnterrs(ROIindx, :) = 100*([FittingResults(ROIbools).T2fit] - ROIinfo(ROIindx).T2)/ROIinfo(ROIindx).T2;
    scatter(ROIindx*ones(Nimage, 1), T2prcnterrs(ROIindx, :), '*')
    hold on
end

boxplot(transpose(T2prcnterrs))
xticks(linspace(1,NROI,NROI))
xticklabels(ROInames)
ylabel('Estimated T2 percentage error')
ylim([-15,15])
grid on
title(opts.disttype)


% === Figure 3: ADC estimation

% Box plot for ADC estimation in each ROI

fig3 = figure;

ROInames  = [];
ADCprcnterrs = zeros(NROI, Nimage);

for ROIindx = 1:NROI

    % Get ROI settings
    ROIname = ROIinfo(ROIindx).ROIname;
    ROInames = [ROInames string(ROIname)];
    ROInum = ROIinfo(ROIindx).ROInum;

    % ROI bools
    ROIbools = ([FittingResults.ROIindx] == ROIindx);

    % Extract T2 measurements error presentages
    ADCprcnterrs(ROIindx, :) = 100*([FittingResults(ROIbools).ADCfit] - ROIinfo(ROIindx).ADC)/ROIinfo(ROIindx).ADC;
    scatter(ROIindx*ones(Nimage, 1), ADCprcnterrs(ROIindx, :), '*')
    hold on
end

boxplot(transpose(ADCprcnterrs))
xticks(linspace(1,NROI,NROI))
xticklabels(ROInames)
ylabel('Estimated ADC percentage error')
ylim([-15,15])
grid on
title(opts.disttype)





% === Figure 4: resnorm

% Box plot for fitting resnrom

fig4 = figure;

ROInames  = [];
resnorms = zeros(NROI, Nimage);

for ROIindx = 1:NROI

    % Get ROI settings
    ROIname = ROIinfo(ROIindx).ROIname;
    ROInames = [ROInames string(ROIname)];
    ROInum = ROIinfo(ROIindx).ROInum;

    % ROI bools
    ROIbools = ([FittingResults.ROIindx] == ROIindx);

    % Extract T2 measurements error presentages
    resnorms(ROIindx, :) = [FittingResults(ROIbools).resnorm];
    scatter(ROIindx*ones(Nimage, 1), resnorms(ROIindx, :), '*')
    hold on
end

boxplot(transpose(resnorms))
xticks(linspace(1,NROI,NROI))
xticklabels(ROInames)
ylabel('Fitting Residual Error')
ylim([0,10e-3])
grid on
title(opts.disttype)


end