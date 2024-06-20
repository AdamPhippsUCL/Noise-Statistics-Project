% MATLAB script to create VERDICT training data from array of T2 and sigma0
% values

% % Define T2 values
% T2min = 50;
% T2max = 100;
% NT2 = 3;
T2s = [50,100,150,200];%linspace(T2min, T2max, NT2);

% % Define sigma0 values
% sigma0min = 0.02;
% sigma0max = 0.05;
% Nsigma0 = 3;
sigma0s = [0.02,0.04,0.06,0.08];%linspace(sigma0min, sigma0max, Nsigma0);


% Define training data folder
TrainDataFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Adaptive Model Fitting\Training Data");

% Define model type
modeltype = 'No VASC VERDICT';

% Define schemename
savescheme = true;
schemename = 'Short Scheme v1';
schemesfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Adaptive Model Fitting\Schemes");

% Define number of training samples
Ntrain = 10000;

savedata = true;

noisetype = 'Ratio';



%% Create data grid
for T2 = T2s
    for sigma0 = sigma0s

        disp(T2)
        disp(sigma0)

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

        pause(1)
       
    end
end
