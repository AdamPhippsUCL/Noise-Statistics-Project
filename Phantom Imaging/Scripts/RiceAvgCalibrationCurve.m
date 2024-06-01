% MATLAB Script to generate calibration curves for averaging Rice
% distributions at different SNR values


% Save folder
folder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Distributions\Calibration Curves";

% Number of samples
Nsample = 25000;

% Rav values
Ravs = [3, 6, 12, 18, 24];

% Define SNR
SNRs = linspace(0.1,10,100);

save([char(folder) '/SNRs.mat'], 'SNRs')



for SNR = SNRs

    disp(SNR)
    v = 1;
    sigma = 1/SNR;
    
    
    %% Define Rice distibution
    
    % Make Rice distribution
    ricedist = makedist('Rician','s',v,'sigma',sigma);
    
    maxsignal = 10;
    signals = linspace(0,maxsignal,501);
    
    pdf = ricedist.pdf(signals);
    
    
    FitVs = zeros(size(Ravs));
    FitSigmas = zeros(size(Ravs));
    
    
    %% Analysis for each Rav value
    
    for RavIndx = 1:length(Ravs)
    
        Rav = Ravs(RavIndx);
        % disp(Rav)
        %% Take samples
        
        % Sample distributions
        samples = zeros(Nsample, 1);
        for indx = 1:Nsample
        
            % Average some signal measurements
            sample = 0;
            for AVindx = 1:Rav
                sample = sample + sampleDistribution(pdf, signals);
            end
            sample = sample/Rav;
        
        
            samples(indx, 1) = sample;
        
        end
        
        
        
        % %% Display original PDF and histogram of samples

        % % Define bin edges
        % binedges = signals;
        % binwidth = binedges(2)-binedges(1);
        % bincenters = binedges(1:end-1)+binwidth/2;
        % 
        % figure;
        % plot(signals, Nsample*binwidth*pdf);
        % hold on
        % H = histogram(samples, binedges);
        % counts = H.Values;
        
        
        
        %% Fit Rician distribution to samples 
        x = (samples) + eps;
        pd = fitdist(x,'Rician');
        
        fitv = pd.s;
        FitVs(RavIndx) = fitv;
    
        fitsigma = pd.sigma;
        FitSigmas(RavIndx) = fitsigma;


        % % Display fitted distribution
        % fitricedist = makedist('Rician','s',fitv,'sigma',fitsigma);
        % fitpdf = fitricedist.pdf(signals);
        % 
        % plot(signals, Nsample*binwidth*fitpdf, Color = 'red')

    
     end
    
    % Make dictionarys
    
    try
        mkdir([char(folder) '/SNR ' num2str(SNR)])
    catch
        disp('')
    end

    MeanDict = dictionary(Ravs, FitVs);
    save([char(folder) '/SNR ' num2str(SNR) '/MeanCalibration.mat'], 'MeanDict' )

    SigmaDict = dictionary(Ravs, FitSigmas/sigma);
    save([char(folder) '/SNR ' num2str(SNR) '/SigmaCalibration.mat'], 'MeanDict' )
    


end



% 
% figure;
% plot(Ravs, FitVs, '*-', color = 'r');
% 
% 
% figure
% plot(Ravs, FitSigmas/sigma, '*-', color = 'b')  





