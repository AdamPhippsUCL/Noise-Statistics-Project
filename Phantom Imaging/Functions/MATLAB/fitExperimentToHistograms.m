function FittingResults = fitExperimentToHistograms(date, experiment, opts)

arguments
    date % Date of phantom imaging session
    experiment % experiment number (details below)

    % OPTIONS
    opts.ROInums = [1, 2, 3, 4, 5, 6, 7, 8] % Specify list of ROI numbers if not all analysed
    opts.disttype = 'Ratio'
    opts.binmin=0
    opts.binmax=2
    opts.nbin=100
    opts.paramtolerance = 0.1 % Fraction tolerance in fitted T2 and ADC values
    opts.sigma0max = 0.25
    opts.DiffDirec = 0 % Diffusion direcion
    opts.imagefolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Imaging Data"
    opts.outputfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Outputs"
end


%% Define folders


% Define main image folder
mainimagefolder = [char(opts.imagefolder) '/'  num2str(date)  '\Images'];
    
% Define numpy image folder
numpyimagefolder = [char(opts.outputfolder) '/'  num2str(date)  '\Images\numpy'];

% Define ROI folder
ROIsfolder = [char(opts.outputfolder) '/'  num2str(date)  '\ROIs'];


% Define ROI folder
ROIfolder = [char(opts.outputfolder) '/'  num2str(date)  '\ROIs'];

% Define ROI values folder
ROIvalsfolder = [char(opts.outputfolder) '/'  num2str(date)  '\ROI values'];

%% Image and ROI details

ImageSettings = createImageSettingStruct(date, experiment);
Nimage = length(ImageSettings);

ROIinfo= createROIinfoStruct(date, ROIfolder);

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

    
    % Check if normalised images are saved, if not: go and save image!
    if ~exist([numpyimagefolder '/' SeriesDescription ], 'dir')
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
            disp(['ROI: ' ROIname ' not saved for image: ' SeriesDescription])
            continue
        end



        %% Extract ROI values
        % ROIvalsfolder = [char(opts.mainfolder) num2str(date)  '\ROI values'];
        
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
            imagefolder = numpyimagefolder,...
            ROIfolder = ROIfolder,...
            savefolder = ROIvalsfolder);

        % Define diffusion orientation
        DiffDirec = opts.DiffDirec;
        
        % Load ROI values
        load([ROIvalsfolder '/' SeriesDescription '/mat/' ROIname '/bval' num2str(bval) '.0_DiffDirec' num2str(DiffDirec) '.mat'])


        %% Make histogram
        
        % Define bin edges and centers
        binmin = opts.binmin;
        binmax = opts.binmax;

        if max(ROIvals)<1
            binmax = 1;
        end
        
        nbin = opts.nbin;
        binedges = linspace(binmin, binmax, nbin+1);
        bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
        binspacing = bincenters(2)-bincenters(1);
        
        f = figure;%('visible','off');
        H = histogram(ROIvals, binedges, FaceColor = [0 0.4470 0.741], FaceAlpha = 0.25, DisplayName = 'ROI values');
        hold on;
        counts = H.Values;
        

        while max(counts) > 0.1*sum(counts)

            nbin = 2*nbin;

            binedges = linspace(binmin, binmax, nbin+1);
            bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
            binspacing = bincenters(2)-bincenters(1);   
            close(f);
            f=figure;
            hold on
            H = histogram(ROIvals, binedges, FaceColor = [0 0.4470 0.741], FaceAlpha = 0.25);
            counts = H.Values;
        end



        % === Fit to histogram
 
        
        % Initial guess
        sigma0guess = std(ROIvals)*(sqrt(2)*exp(TE/T2));
        fdguess = mean(ROIvals);

        beta0guess = [0.018, T2, exp(-ADC*bval)];

        % Bounds
        % ADClb = ROIinfo(ROIindx).ADClb;
        % ADCub = ROIinfo(ROIindx).ADCub;
        % T2lb = ROIinfo(ROIindx).T2lb;
        % T2ub = ROIinfo(ROIindx).T2ub;

        ADClb = ADC - 0.5e-4;
        ADCub = ADC + 0.5e-4;

        T2lb = T2-100;
        T2ub = T2+100;
        
        if strcmp(opts.disttype, 'Rice')
            lb = [0.001, T2-eps, exp(-(ADCub)*bval)];
            ub = [opts.sigma0max, T2+eps, exp(-(ADClb)*bval)];
        elseif strcmp(opts.disttype, 'Ratio')
            lb = [0.001, T2lb, exp(-(ADCub)*bval)];
            ub = [opts.sigma0max, T2ub, exp(-(ADClb)*bval)];
        end

        % if 0.999>b*exp(-ADC*bval)
        %     ub = [opts.sigma0max, b*T2, b*exp(-ADC*bval)];
        % else
        %     ub = [opts.sigma0max, b*T2,0.9999];
        % end


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
        
        xlabel('Normalized signal')
        ylabel('Counts')
        if strcmp(opts.disttype, 'Ratio')
            plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = [opts.disttype ' distribution'], color = "#D95319");
        else
            plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = [opts.disttype ' distribution'], color = "#7E2F8E");
        end

        ax = gca();
        ax.FontSize = 12;
        legend;
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

colordict = dictionary([1,2,3,4,5,6,7,8], [	"#0072BD", 	"#D95319",	"#EDB120", 	"#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#e981cd"]);
concentrations = ["50% PVP", "40% PVP", "30% PVP", "20% PVP", "10% PVP", "5% PVP", "2.5% PVP", "1mM NiCl_2"];

switch experiment

    case 1
        
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
            scatter(ROIindx*ones(Nimage, 1), Sigma0s(ROIindx, :), '*', MarkerEdgeColor = colordict(ROIindx), DisplayName = concentrations(ROIindx)) 
            hold on
        
        end
        
        boxplot(transpose(Sigma0s))
        xticks(linspace(1,NROI,NROI))
        xticklabels(ROInames)
        ylabel('Estimated \sigma_0')
        % xlabel('Concentration')
        ylim([0, 0.2]);%opts.sigma0max])
        grid on
        title(opts.disttype)
        legend;
        ax = gca();
        ax.FontSize = 12;
        
        % === Figure 2: T2 estimation
        
        % Box plot for ADC estimation in each ROI
        
        fig2 = figure;
        
        ROInames  = [];
        T2s = zeros(NROI, Nimage);
        T2errs = zeros(NROI, Nimage);
        
        for ROIindx = 1:NROI
        
            % Get ROI settings
            ROIname = ROIinfo(ROIindx).ROIname;
            ROInames = [ROInames string(ROIname)];
            ROInum = ROIinfo(ROIindx).ROInum;
        
            % ROI bools
            ROIbools = ([FittingResults.ROIindx] == ROIindx);
        
            % % Extract T2 measurements error presentages
            % T2errs(ROIindx, :) = ([FittingResults(ROIbools).T2fit] - ROIinfo(ROIindx).T2);%/ROIinfo(ROIindx).T2;
            % scatter(ROIindx*ones(Nimage, 1), T2errs(ROIindx, :), '*')
            % hold on

            trueT2 = ROIinfo(ROIindx).T2;
            T2s(ROIindx, : ) = ([FittingResults(ROIbools).T2fit]);
            scatter(ROIindx*ones(Nimage, 1), T2s(ROIindx, :), '*', MarkerEdgeColor = colordict(ROIindx), DisplayName = concentrations(ROIindx)) 
            hold on
            scatter([ROIindx+0.4], [trueT2], MarkerEdgeColor = 'black', MarkerFaceColor='black', Marker = 'o', HandleVisibility = 'off')

        end
        
        boxplot(transpose(T2s))
        xticks(linspace(1,NROI,NROI))
        xticklabels(ROInames)
        ylabel('Estimated T2')
        % xlabel('Concentration')
        ylim([0,1000])
        grid on
        title(opts.disttype)
        % legend;
        ax = gca();
        ax.FontSize = 12;
        
        % === Figure 3: ADC estimation
        
        % Box plot for ADC estimation in each ROI
        
        fig3 = figure;
        
        ROInames  = [];
        ADCs = zeros(NROI, Nimage);
        ADCerrs = zeros(NROI, Nimage);
        fdprcnterrs = zeros(NROI, Nimage);
        
        for ROIindx = 1:NROI
        
            % Get ROI settings
            ROIname = ROIinfo(ROIindx).ROIname;
            ROInames = [ROInames string(ROIname)];
            ROInum = ROIinfo(ROIindx).ROInum;
        
            % ROI bools
            ROIbools = ([FittingResults.ROIindx] == ROIindx);
        
            % % Extract ADC measurements error presentages
            % ADCprcnterrs(ROIindx, :) = 100*([FittingResults(ROIbools).ADCfit] - ROIinfo(ROIindx).ADC)/ROIinfo(ROIindx).ADC;
            % scatter(ROIindx*ones(Nimage, 1), ADCprcnterrs(ROIindx, :), '*')
            % hold on

            % % Extract fd measurements error precentages
            % theseFittingResults = FittingResults(ROIbools);
            % ADCfits = [theseFittingResults.ADCfit];
            % 
            % trueADC = ROIinfo(ROIindx).ADC;
            % 
            % for imgindx = 1:Nimage
            %     thisADC = ADCfits(imgindx);
            %     thisimgsettings = theseFittingResults(imgindx).ImageSettings;
            %     thisbval = thisimgsettings.bval;
            % 
            % 
            %     thisfd = exp(-thisADC*thisbval);
            %     truefd = exp(-trueADC*thisbval);
            % 
            %     fdprcnterr = 100*(thisfd-truefd)/truefd;
            % 
            %     fdprcnterrs(ROIindx, imgindx) = fdprcnterr;
            % 
            % 
            % end
            
            % ADCerrs(ROIindx, :) = ([FittingResults(ROIbools).ADCfit] - ROIinfo(ROIindx).ADC);%/ROIinfo(ROIindx).ADC;
            % scatter(ROIindx*ones(Nimage, 1), ADCerrs(ROIindx, :), '*')
            % hold on

            trueADC = ROIinfo(ROIindx).ADC;
            ADCs(ROIindx, :) = ([FittingResults(ROIbools).ADCfit]);
            scatter(ROIindx*ones(Nimage, 1), ADCs(ROIindx, :), '*', MarkerEdgeColor = colordict(ROIindx), DisplayName = concentrations(ROIindx))
            hold on
            scatter([ROIindx+0.4], [trueADC], MarkerEdgeColor = 'black', MarkerFaceColor='black', Marker = 'o', HandleVisibility='off')

        end
        
        % boxplot(transpose(fdprcnterrs))
        boxplot(transpose(ADCs))
        xticks(linspace(1,NROI,NROI))
        xticklabels(ROInames)
        ylabel('Estimated D (mm^2/s)')
        % xlabel('Concentration')
        ylim([0,1.2e-3])
        grid on
        title(opts.disttype)
        % legend;
        ax = gca();
        ax.FontSize = 12;
        
        
        
        
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
        % ylim([0,10e-3])
        grid on
        title(opts.disttype)
        ax = gca();
        ax.FontSize = 12;
        


    case 2

        figure;

        % Scatter plot for each ROI
        for ROIindx = opts.ROInums

            bools = [FittingResults.ROIindx] == ROIindx;
            sigma0s = [FittingResults(bools).sigma0];
            Ravs = [ImageSettings(:).Nav_ratio];

            plot(Ravs, sigma0s, '*-', DisplayName = ['ROI ' num2str(ROIindx)], Color = colordict(ROIindx))
            hold on
            xlabel('N_{b>0}')
            ylabel('\sigma_0')
            ylim([0,0.1])
            legend(NumColumns = 2)
            ax = gca();
            ax.FontSize = 12;
        
        end

        figure;

        % Scatter plot for each ROI
        for ROIindx = opts.ROInums

            trueT2 = ROIinfo(ROIindx).T2;
            bools = [FittingResults.ROIindx] == ROIindx;
            T2s = [FittingResults(bools).T2fit];
            Ravs = [ImageSettings(:).Nav_ratio];

            plot(Ravs, T2s, '*-', DisplayName = ['ROI ' num2str(ROIindx)], Color = colordict(ROIindx))
            hold on
            plot(Ravs, trueT2*ones(size(Ravs)), '--', Color = colordict(ROIindx))
            xlabel('N_{b>0}')
            ylabel('T2 (ms)')
            ylim([0,1000])
            xlim([1, 25])
            ax = gca();
            ax.FontSize = 12;       
            % legend
        end


        figure;

        % Scatter plot for each ROI
        for ROIindx = opts.ROInums

            trueADC = ROIinfo(ROIindx).ADC;
            bools = [FittingResults.ROIindx] == ROIindx;
            ADCs = [FittingResults(bools).ADCfit];
            Ravs = [ImageSettings(:).Nav_ratio];

            plot(Ravs, ADCs, '*-', DisplayName = ['ROI ' num2str(ROIindx)], Color = colordict(ROIindx))
            hold on
            plot(Ravs, trueADC*ones(size(Ravs)), '--', Color = colordict(ROIindx))
            % plot(Ravs, )
            xlabel('N_{b>0}')
            ylabel('D (mm^2/s)')
            ylim([0,1.2e-3])
            xlim([1, 25])
            ax = gca();
            ax.FontSize = 12;       
            % legend
        end
       
    case 3

        figure;

        % Scatter plot for each ROI
        for ROIindx = opts.ROInums

            bools = [FittingResults.ROIindx] == ROIindx;
            sigma0s = [FittingResults(bools).sigma0];
            NSAs = [ImageSettings(:).NSA];

            plot(1./sqrt(NSAs), (sigma0s-sigma0s(end))/(sigma0s(1)-sigma0s(end)), '*-', DisplayName = ['ROI ' num2str(ROIindx)])
            hold on
            xlabel('NSA^{-1/2}')
            ylabel('\sigma_0')
            % ylim([0,0.1])
            legend(NumColumns = 2)

        end

end

