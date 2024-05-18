% MATLAB function to evaluate fIC fitting performance for different T2 and
% sigma0 combinations
function [FittingBiases, FittingVariances] = NS_EvaluateFitting(opts)

arguments

    % Noise parameters 
    opts.T2s = [40, 80, 120, 40, 80, 120, 40, 80, 120]
    opts.sigma0s = [0.01, 0.01, 0.01, 0.0275, 0.0275, 0.0275, 0.045, 0.045, 0.045]

    % parameters of constant ratio distribution
    opts.T2const = 10000;
    opts.sigma0const = 0.02;

    % Simulation
    opts.modeltype = 'RDI v1.3'
    opts.schemename = 'Short Scheme v1'
    opts.fittingtype = 'MLP'
    opts.schemesfolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\General Code\Model Fitting\MLP\My Tests\Schemes'
    opts.fICs = [0.1, 0.9] % min and max fIC values to use in simulation
    opts.fVASCs = [0, 0.1] % min and max fVASC values to use in simulation
    opts.Rs = linspace(0.1, 15.1, 17) % radii values to use in simulation
    opts.randtype = 'normal' % How fR distibution is randomised (uniform or Gaussian distribution)
    opts.muRs = [6,9] % mu and sigma ranges of radii distribution 
    opts.sigmaRs = [1,3]

    % Fitting
    opts.Niter = 100
    opts.Nrep = 100

    opts.figpath= 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Figures'


end

% Number of experiments
Nprotocols = length(opts.T2s);


% Initialise array for results
FittingBiases = zeros(2, Nprotocols, opts.Niter);
FittingVariances = zeros(2, Nprotocols, opts.Niter);

% Summary stats
stats = struct();


%% Evaluate fitting performance of each protocol

% Track simulated parameters
simfICs = zeros(opts.Niter,1);
simfRs = zeros(opts.Niter,length(opts.Rs));

% == Iterate over randomised tissue parameters

for iterIndx = 1:opts.Niter

    iterIndx
    
    % == Find randomised paramater set (ALWAYS SIMULATING WITH ORIGINAL
    % VERDICT MODEL)

    % fIC, fEES, fVASC
    fIC = opts.fICs(1) + ( opts.fICs(2) - opts.fICs(1) )*rand();

    if 1-fIC > opts.fVASCs(2)
        fVASC = opts.fVASCs(1) + (opts.fVASCs(2)-opts.fVASCs(1))*rand();
    else
        fVASC = opts.fVASCs(1) + (1-fIC-opts.fVASCs(1))*rand();
    end

    fEES = 1-fIC-fVASC;

    switch opts.randtype
    
        case 'uniform'

            % Radii distribution (FROM UNIFORM CURRENTLY)
            fRs = rand(length(opts.Rs), 1);

        case 'normal'
        
            % Radii distribution (Normal)
            muR = opts.muRs(1) + ( opts.muRs(2) - opts.muRs(1) )*rand();
            sigmaR = opts.sigmaRs(1) + ( opts.sigmaRs(2) - opts.sigmaRs(1) )*rand();
            fRs = normpdf(opts.Rs, muR, sigmaR);

    end

    % Scale by fIC
    fRs = (fIC/sum(fRs))*fRs;

    simfICs(iterIndx) = fIC;
    simfRs(iterIndx,:) = fRs;


    % == Evaluate fitting at this tissue parameter set for each protocol

    for protIndx = 1:Nprotocols

        % modeltype = modeltypes{protIndx};
        % schemename = schemenames{protIndx};
        % fittingtype = fittingtypes{protIndx};
        modeltype = opts.modeltype;
        schemename = opts.schemename;
        fittingtype = opts.fittingtype;

        % Load scheme
        load([opts.schemesfolder '/' schemename '.mat']);

      
        % Evaluate fitting (CONSTANT RATIO DISTRIBUTION)
        [bias, variance] = nsevaluatefitting( ...
            scheme, ...
            schemename,...
            modeltype, ...
            fittingtype,...
            fIC = fIC,...
            fEES = fEES,...
            Rs = opts.Rs,...
            fRs = fRs,...
            T2 = opts.T2s(protIndx),...
            NoiseSigma = opts.sigma0s(protIndx),...
            Nrep = opts.Nrep,...
            sigma0train = opts.sigma0const,...
            T2train = opts.T2const...
            );

        FittingBiases(1, protIndx, iterIndx) = bias;
        FittingVariances(1, protIndx, iterIndx) = variance;

        % Evaluate fitting (NOISE ADAPTIVE RATIO DISTRIBUTION)
        [bias, variance] = nsevaluatefitting( ...
            scheme, ...
            schemename,...
            modeltype, ...
            fittingtype,...
            fIC = fIC,...
            fEES = fEES,...
            Rs = opts.Rs,...
            fRs = fRs,...
            T2 = opts.T2s(protIndx),...
            NoiseSigma = opts.sigma0s(protIndx),...
            Nrep = opts.Nrep,...
            sigma0train = opts.sigma0s(protIndx),...
            T2train = opts.T2s(protIndx)...
            );
    
    
        FittingBiases(2, protIndx, iterIndx) = bias;
        FittingVariances(2, protIndx, iterIndx) = variance;

    end

end




n=size(FittingVariances,3);
nprot = size(FittingBiases,2);
maxbias = max(FittingBiases(:));
minbias = min(FittingBiases(:));
maxvar = max(FittingVariances(:));

for protindx = 1:nprot

    fig1 = figure;
    scatter(ones(n,1), squeeze(FittingBiases(1,protindx,:)), '*', DisplayName = 'Rice', MarkerEdgeAlpha=0.3)
    hold on
    scatter(2*ones(n,1), squeeze(FittingBiases(2,protindx,:)), '*', DisplayName = 'Ratio', MarkerEdgeAlpha=0.3)    
    boxplot(transpose(squeeze(FittingBiases(:,protindx,:)) ))
    ylim([1.1*minbias, 1.1*maxbias])
    ylabel('f_{IC} fitting bias')
    title( ['T2: ' num2str(opts.T2s(protindx)) ', \sigma_0: ' num2str(opts.sigma0s(protindx))])
    legend
    saveas(fig1,  [opts.figpath '/Bias T2 ' num2str(opts.T2s(protindx)) ', sigma0 ' num2str(opts.sigma0s(protindx)) '.png'])

    fig2 = figure;
    scatter(ones(n,1), squeeze(FittingVariances(1,protindx,:)), '*', MarkerEdgeAlpha=0.3)
    hold on
    scatter(2*ones(n,1), squeeze(FittingVariances(2,protindx,:)), '*', MarkerEdgeAlpha=0.3)    
    boxplot(transpose(squeeze(FittingVariances(:,protindx,:)) ))
    ylim([0, 1.1*maxvar])
    ylabel('f_{IC} fitting variance')
    title( ['T2: ' num2str(opts.T2s(protindx)) ', \sigma_0: ' num2str(opts.sigma0s(protindx))])
    saveas(fig2,  [opts.figpath '/Variance T2 ' num2str(opts.T2s(protindx)) ', sigma0 ' num2str(opts.sigma0s(protindx)) '.png'] )
end


end