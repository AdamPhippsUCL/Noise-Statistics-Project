% MATLAB script to train VERDICT MLPs from array of T2 and sigma0
% values


% % Define T2 values
% T2min = 50;
% T2max = 125;
% NT2 = 4;
T2s = [50,100,150,200];%40,60,80,100];%linspace(T2min, T2max, NT2);

% % Define sigma0 values
% sigma0min = 0.02;
% sigma0max = 0.06;
% Nsigma0 = 3;
sigma0s = [0.02,0.04,0.06,0.08];%, 0.033, 0.05];%linspace(sigma0min, sigma0max, Nsigma0);

noisetype = 'Ratio';

% Define training data folder
TrainDataFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Adaptive Model Fitting\Training Data");

% Model folder
ModelFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Adaptive Model Fitting\MLP Models");

% Python folder
pythonfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\General VERDICT Code\General-VERDICT-Code\MLP\Python");

% Define model type
modeltype = 'No VASC VERDICT';

% Define schemename
schemename = 'Short Scheme v1';

schemesfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\Adaptive Model Fitting\Schemes");


% Save sigma0 and T2 values
save([ModelFolder '/' modeltype '/' schemename '/' noisetype '/T2s.mat'], "T2s")
save([ModelFolder '/' modeltype '/' schemename '/' noisetype '/sigma0s.mat'], "sigma0s")


%% Train MLP grid
for T2 = T2s
    for sigma0 = sigma0s

        % update model and training data folders
        % this_schemename = [schemename '/T2_' num2str(T2) '/sigma_' num2str(sigma0)];
        trainMLP( ...
            modeltype, ...
            schemename, ...
            noisetype=noisetype, ...
            T2train = T2, ...
            sigma0train = sigma0, ...
            TrainDataFolder=TrainDataFolder, ...
            ModelFolder=ModelFolder,...
            pythonfolder = pythonfolder);
    end
end