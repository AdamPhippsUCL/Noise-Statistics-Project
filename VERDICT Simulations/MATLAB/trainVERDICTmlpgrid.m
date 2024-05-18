% MATLAB script to train VERDICT MLPs from array of T2 and sigma0
% values


% Define T2 values
T2min = 50;
T2max = 125;
NT2 = 4;
T2s = linspace(T2min, T2max, NT2);

% Define sigma0 values
sigma0min = 0.02;
sigma0max = 0.06;
Nsigma0 = 3;
sigma0s = [0.067];%, 0.033, 0.05];%linspace(sigma0min, sigma0max, Nsigma0);



% Define training data folder
TrainDataFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Training Data");

% Model folder
ModelFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\MLP Models");


% Define model type
modeltype = 'RDI v1.3';

% Define schemename
schemename = 'Short Scheme NS';
schemesfolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Schemes';



%% Train MLP grid
for T2 = T2s
    for sigma0 = sigma0s


        % update model and training data folders
        this_schemename = [schemename '/T2_' num2str(T2) '/sigma_' num2str(sigma0)];
        trainMLP(modeltype, this_schemename, TrainDataFolder=TrainDataFolder, ModelFolder=ModelFolder)
    end
end