% MATLAB function to generate training data 

function [params, Signals] = createVERDICTdataNew(modeltype, schemename, opts)

arguments

    modeltype % specify model type
    schemename % specify scheme name

    % == OPTIONS

    opts.Ntrain = 100% Number of training samples

    
    % Noise
    opts.noisetype = 'Ratio'
    opts.sigma0 = 0.05 % 
    opts.T2 = 100 % (ms) 
    opts.TEconst = 22 % (ms) TE = delta + Delta + TEconst

    % Parameter ranges (restricted to better represent real tissue)
    opts.fICs = [0.1, 0.9]
    opts.fVASCs = [0, 0.2]

    % == VERDICT model
    opts.Rs = linspace(0.1,15.1,17)
    opts.randtype = 'normal'
    opts.randmuRs = [5,10] % range of fR distribution means
    opts.randsigmaRs = [1,3] % range of fR distribution sigmas

    % == RDI model
    opts.muRs = [5, 7.5, 10]
    opts.sigmaRs = [2,2,2]
    opts.muRrange = [5,10] % Range of R/muR for RDI v1.3 and v1.4
    opts.sigmaRrange = [1,3] % Range of sigmaR for RDI v1.4



    % Save scheme
    opts.savescheme = true;
    opts.schemeParentFolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\General Code\Model Fitting\MLP\My Tests\Schemes'

    % Save data
    opts.savedata = false
    opts.dataParentFolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\General Code\Model Fitting\MLP\My Tests\Training Data'

end

% STEPS

% 1. Build specified scheme
% 2. Create set of randomised model parameters (for specified model)
% 3. Simulate signal over scheme for each parameter set (distribution sampling)
% 4. Save training data

% =========================================================================

%% 1. Build scheme

% [delta, DELTA, b, TE, TR, NSA, Rav]  
% (NSA = number of signal averages for b=0 image)
% (Rav = additional averaging factor for diffusion weighted image)

switch schemename

 
    case 'Short Scheme v1'

        V01 = [1, 2, 0, 76, 4000, 4, 1];
        V1  = [20, 41, 1800, 76, 4000, 4, 3*3]; % 22.06.24 correct (Rav = DWI*AvgHighB)
        V02 = [1, 2, 0, 62, 4000, 4, 1];
        V2 = [12, 35, 1000, 62, 4000, 4, 3*3];

        Vs = [...
            V01; V1;...
            V02; V2;...
            ];

        % build scheme
        scheme = BuildSchemeNew(Vs);

end


if opts.savescheme
    save([opts.schemeParentFolder '\' schemename '.mat'], 'scheme')
end


%% 2. Randomise model paramaters

% == Specify number of parameters
switch modeltype

    case 'Original VERDICT'

        % Define radii
        Rs = opts.Rs;
        nR = length(Rs);

        % Define ncompart
        ncompart = 2;

        % Number of parameters
        nparam = nR + ncompart;

    case 'No VASC VERDICT'

        % Define radii
        Rs = opts.Rs;
        nR = length(Rs);

        % Define ncompart
        ncompart = 1;

        % Number of parameters
        nparam = nR + ncompart;


    case 'RDI v1.3'

        ncompart = 1;
        nparam = 2 + ncompart;


        

end


% == Randomise parameters

% Initialise parameter array
params = zeros(opts.Ntrain, nparam); 

% Loop to fill array
for paramIndx = 1:opts.Ntrain

    % == Intracellular

%     % Randomise fIC
%     fIC = opts.fICs(1) + (opts.fICs(2)-opts.fICs(1))*rand();

    %% == INVESTIGATING FIC randomisation distribution!
    fICs = linspace(0,1,100);
    % pdf = ones(size(fICs));
    % pdf = 2*abs(1-2*fICs);
    pdf = (2/3)*(1 + 1*(cos(pi*fICs).^2));
%     pdf = (2/3)*(1 + sin(pi*fICs).^2);
    fIC = sampleDistribution(pdf, fICs);

    % number of IC compartments (induvidual radii or distributions)
    nIC = nparam - ncompart;


    switch modeltype
        
        case 'Original VERDICT'

            switch opts.randtype

                % fRs randomly generated from uniform distribution over Rs
                case 'uniform'

                    % uniform distribution
                    fs = rand(1,nIC);
                
                    % Normalise
                    fs = fIC*fs/sum(fs);
    
                
                % fRs generates as random uniform distribution
                case 'normal'
                    
                    randmuR = opts.randmuRs(1) + (opts.randmuRs(2) - opts.randmuRs(1))*rand();
                    randsigmaR = opts.randsigmaRs(1) + (opts.randsigmaRs(2) - opts.randsigmaRs(1))*rand();

                    fs = normpdf(Rs, randmuR, randsigmaR );

                    % Normalise
                    fs = fIC*fs/sum(fs);

            end


        case 'No VASC VERDICT'

            switch opts.randtype

                case 'uniform'
                    
                    % uniform distribution
                    fs = rand(1,nIC);
                
                    % Normalise
                    fs = fIC*fs/sum(fs);

                case 'normal'
                    
                    randmuR = opts.randmuRs(1) + (opts.randmuRs(2) - opts.randmuRs(1))*rand();
                    randsigmaR = opts.randsigmaRs(1) + (opts.randsigmaRs(2) - opts.randsigmaRs(1))*rand();

                    fs = normpdf(Rs, randmuR, randsigmaR );

                    % Normalise
                    fs = fIC*fs/sum(fs);

            end       



        case 'RDI v1.3'

            fIC = rand();
            R = opts.muRrange(1) + (opts.muRrange(2)-opts.muRrange(1))*rand();

            fs = [fIC, R];



    end



    % == Remaining volume fractions

    % ncompart = 2
    if ncompart == 2

        if 1-fIC > opts.fVASCs(2)
            fVASC = opts.fVASCs(1) + (opts.fVASCs(2)-opts.fVASCs(1))*rand();
        else
            fVASC = opts.fVASCs(1) + (1-fIC-opts.fVASCs(1))*rand();
        end

        fEES = 1-fIC-fVASC;

        fend = [fEES, fVASC];
    

    % ncompart = 1    
    else
        fEES = 1-fIC;
        fVASC = 0;
        fend = fEES;
    end

    % == Append volume fractions to array
    params(paramIndx, :) = [fs, fend];


end



%% 3. Simulate signal over parameter sets

% Initialise signal array
Signals = zeros(opts.Ntrain, length(scheme));

switch modeltype

    % case 'Original VERDICT'
    % 
    %     for paramIndx = 1:opts.Ntrain
    % 
    %         paramIndx/opts.Ntrain
    % 
    %         tps = params(paramIndx,:);
    % 
    %         % Simulate signal distributions
    %         signals = simulateSignalsOverScheme(...
    %         scheme,...
    %         modeltype,...
    %         fIC = sum(tps(1:end-2)),...
    %         fEES = tps(end-1),...
    %         fVASC = tps(end),...
    %         fRs = tps(1:end-2),...
    %         Rs = Rs,...
    %         sigma0 = opts.sigma0,...
    %         T2 = opts.T2,...
    %         TEconst = opts.TEconst,...
    %         noisetype=opts.noisetype...
    %         );
    % 
    %         % f=figure;
    %         % scatter([scheme(:).bval], signals)
    %         % close(f);
    % 
    %         Signals(paramIndx,:) = signals;
    % 
    % 
    %      end

    case 'No VASC VERDICT'

        for paramIndx = 1:opts.Ntrain

            paramIndx/opts.Ntrain

            tps = params(paramIndx,:);
    
            % Simulate signal distributions
            signals = simulateSignalsOverSchemeNew(...
            scheme,...
            modeltype,...
            fIC = sum(tps(1:end-1)),...
            fEES = tps(end),...
            fVASC = 0,...
            fRs = tps(1:end-1),...
            Rs = Rs,...
            sigma0 = opts.sigma0,...
            T2 = opts.T2,...
            TEconst = opts.TEconst,...
            noisetype=opts.noisetype...
            );

            Signals(paramIndx,:) = signals;


        end

    % 
    % case 'RDI v1.3'
    % 
    % 
    %     for paramIndx = 1:opts.Ntrain
    % 
    %         % paramIndx/opts.Ntrain
    % 
    %         tps = params(paramIndx,:);
    % 
    %         % Simulate signal distributions 
    %         signals = simulateSignalsOverScheme(...
    %         scheme,...
    %         modeltype,...
    %         fIC = tps(1),...
    %         fEES = tps(end),...
    %         fVASC = 0,...
    %         muR = tps(2),...
    %         sigma0 = opts.sigma0,...
    %         T2 = opts.T2,...
    %         TEconst = opts.TEconst,...
    %         noisetype=opts.noisetype...
    %         );
    % 
    %         Signals(paramIndx,:) = signals;
    % 
    %     end



end



%% 4. Save training data
if opts.savedata

    % Define output folder
    outputFolder = join([opts.dataParentFolder '/' modeltype '/' schemename '/' opts.noisetype '/T2_' num2str(opts.T2) '/sigma_' num2str(opts.sigma0)], "");
    
    if ~exist(outputFolder, "dir")
       mkdir(outputFolder)
    end

    save([outputFolder '/params.mat'], 'params');
    save([outputFolder '/signals.mat'], 'Signals');

    % Save META data structure
    Meta = struct();
    Meta.DateTime = datetime();
    Meta.Ntrain = opts.Ntrain;
    Meta.sigma0 = opts.sigma0;
    Meta.T2 = opts.T2;
    Meta.TEconst = opts.TEconst;
    
    Meta.modeltype = modeltype;
    Meta.scheme = scheme;
    Meta.schemename = schemename;
    Meta.Rs = opts.Rs;
    Meta.randtype = opts.randtype;
    Meta.fICpdf = pdf;
    Meta.randmuRs = opts.randmuRs;
    Meta.randsigmaRs = opts.randsigmaRs; 
    Meta.muRs = opts.muRs;
    Meta.sigmaRs = opts.sigmaRs;
    Meta.noisetype = opts.noisetype;

    save([outputFolder '/Meta.mat'], 'Meta');

end

