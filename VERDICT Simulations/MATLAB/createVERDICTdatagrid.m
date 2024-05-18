% MATLAB script to create VERDICT training data from array of T2 and sigma0
% values

% Define T2 values
T2min = 50;
T2max = 100;
NT2 = 3;
T2s = linspace(T2min, T2max, NT2);

% Define sigma0 values
sigma0min = 0.02;
sigma0max = 0.05;
Nsigma0 = 3;
sigma0s = [0.067, 0.033, 0.05];%linspace(sigma0min, sigma0max, Nsigma0);


% Define training data folder
TrainDataFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Training Data");

% Define model type
modeltype = 'RDI v1.3';

% Define schemename
savescheme = true;
schemename = 'NS Simulation Scheme';
schemesfolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Schemes';

% Define number of training samples
Ntrain = 10000;

savedata = true;

noisetype = 'Ratio';

%% Create data grid
for T2 = T2s
    for sigma0 = sigma0s

        createVERDICTdata( ...
            modeltype, ...
            schemename, ...
            T2 = T2,...
            sigma0 = sigma0,...
            Ntrain=Ntrain,...
            savescheme=savescheme,...
            schemeParentFolder=schemesfolder,...
            savedata = true,...
            dataParentFolder = TrainDataFolder,...
            noisetype=noisetype);
       
    end
end
