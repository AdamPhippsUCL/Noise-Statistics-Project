% Matlab script to make simulated image with different T2 and sigma0 values


% sigma0
sigma0 = 0.067;

% T2 values
T2min = 50;
T2max = 100;
NT2 = 3;
T2s = linspace(T2min, T2max, NT2);

% fIC
fICmin = 0.1;
fICmax = 0.9;
NfIC = 9;
fICs = linspace(fICmin, fICmax, NfIC);

% Radii 
Rmin = 6;
Rmax = 9;

%% Make image

% grid size
ngrid = 25;

% T2 grid
T2grid = zeros(NT2*ngrid, NfIC*ngrid);

for T2indx = 1:NT2
    for gridindx = 1:ngrid
        gridindx
        T2grid( (T2indx-1)*ngrid + gridindx, :) = T2s(T2indx);
    end
end

% fIC grid
fICgrid = zeros(NT2*ngrid, NfIC*ngrid);

for fICindx = 1:NfIC
    for gridindx = 1:ngrid
        fICgrid( :, (fICindx-1)*ngrid + gridindx) = fICs(fICindx);
    end
end

% Display simulated fICs
fig = figure;
imshow(fICgrid, [0,1], 'InitialMagnification', 100);
c = colorbar;
ylabel(c, 'f_{IC}')
h = gca;
h.Visible = 'On';
yticks([round(ngrid/2) + ngrid*(0:NT2-1) ])
yticklabels(T2s)
ylabel('T2 (ms)')
xticks([round(ngrid/2) + ngrid*(0:NfIC-1) ])
xticklabels(fICs)
xlabel('f_{IC}')
title('Simulated fIC')
saveas(fig, 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Figures\simfIC.fig')

% Radii (randomise)
Rgrid = Rmin + rand(NT2*ngrid, NfIC*ngrid)*(Rmin - Rmax);



%% Simulate VERDICT signals over image

modeltype = 'RDI v1.3';
schemename = 'Short Scheme NS';
schemesfolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Schemes';
load([schemesfolder '/' schemename '.mat'])
nscheme = length(scheme);

% Signals 
Signals = zeros(NT2*ngrid, NfIC*ngrid, nscheme);

% Simulate signals
for jndx = 1:size(T2grid,1)
    jndx
    for indx = 1:size(fICgrid, 2)
        
        T2 = T2grid(jndx, indx);
        fIC = fICgrid(jndx, indx);
        fEES = 1-fIC;
        R = Rgrid(jndx, indx);

        signals = simulateSignalsOverScheme( ...
            scheme, ...
            modeltype,...
            fIC = fIC,...
            fEES=fEES,...
            muR = R,...
            sigma0 = sigma0,...
            T2 = T2);

        Signals(jndx,indx,:)=signals;



    end
end


%% Fitting

noisetypes = {'Ratio', 'Rice'};
sigma0rice = 0.05;

BiasResults = zeros(length(noisetypes), length(T2s), length(fICs));
VarianceResults = zeros(length(noisetypes),length(T2s), length(fICs));

for nindx = 1:length(noisetypes)

    noisetype = noisetypes{nindx};
    
    switch noisetype
        case 'Ratio'
        % Possible T2 and sigma0 values for Ratio
        possibleT2s = [50, 75, 100, 125];
        possiblesigma0s = [0.02, 0.033, 0.05];
    
        case 'Rice'
        % Possible T2 and sigma0 values for Rice
        possibleT2s = [10000];
        possiblesigma0s = [sigma0rice];
    end


    % Do the fitting ting!
    [fIC, fEES, fVASC, R, rmse] = ns_verdict_MLP_fit( ...
        schemename, ...
        modeltype, ...
        Signals, ...
        noisetype=noisetype,...
        sigma0train = sigma0, ...
        T2train  = T2grid, ...
        possibleT2s=possibleT2s,...
        possiblesigma0s=possiblesigma0s,...
        schemesfolder=schemesfolder...
        );

    % Display
    fig = figure;
    imshow(fIC, [0 1], 'InitialMagnification', 100);
    c = colorbar;
    ylabel(c, 'f_{IC}')
    h = gca;
    h.Visible = 'On';
    yticks([round(ngrid/2) + ngrid*(0:NT2-1) ])
    yticklabels(T2s)
    ylabel('T2 (ms)')
    xticks([round(ngrid/2) + ngrid*(0:NfIC-1) ])
    xticklabels(fICs)
    xlabel('f_{IC}')

    if strcmp(noisetype, 'Rice')
        title(['Simulation: \sigma_0 = ' num2str(sigma0) '; Fitting: Rice \sigma_0 = ' num2str(sigma0)])
        saveas(fig, ['C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Figures\fIC ' noisetype ' (Rice sigma0 ' num2str(sigma0rice) ') sigma0 ' num2str(sigma0) '.fig'])
    else
        title(['Simulation: \sigma_0 = ' num2str(sigma0) '; Fitting: ' char(noisetypes(nindx)) ' \sigma_0 = ' num2str(sigma0)])
        saveas(fig, ['C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Figures\fIC ' noisetype ' sigma0 ' num2str(sigma0) '.fig'])
   
    end


    % Bias and variance as function of fIC
    for T2indx = 1:length(T2s)
        T2val = T2s(T2indx);
        for fICindx = 1:length(fICs)
            fICval = fICs(fICindx);
            mask = and(fICgrid == fICval, T2grid == T2val);
            fICvals = fIC.*mask;
            fICvals = fICvals(fICvals~=0);
            fICvals = fICvals(:);
            bias = mean(fICvals-fICval);
            BiasResults(nindx, T2indx, fICindx) = bias;
            variance = var(fICvals)
            VarianceResults(nindx, T2indx, fICindx) = variance;
        end
    end
end


%% Save Bias and variance results
ResultsFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\MATLAB\Results");

ResultsMeta = struct();
ResultsMeta.noisetypes = noisetypes;
ResultsMeta.T2s = T2s;
ResultsMeta.fICs = fICs;

save([ResultsFolder '/Meta_sigma_' num2str(sigma0) '.mat' ], "ResultsMeta");
save([ResultsFolder '/Bias_sigma_' num2str(sigma0) '.mat' ], "BiasResults");
save([ResultsFolder '/Variance_sigma_' num2str(sigma0) '.mat' ], "VarianceResults");

%% Display bias and variance results

noisemarkers = {'-*','-o', '-p'};
T2colors = {"#0072BD", 	"#77AC30", "#D95319", "#4DBEEE"};

% Bias figure
fig = figure;
title('Bias')
for nindx = 1:length(noisetypes)
    for T2indx = 1:length(T2s)
        plot(fICs, squeeze(BiasResults(nindx,T2indx,:)), noisemarkers{nindx}, color = T2colors{T2indx}, DisplayName=['T2: ' num2str(T2s(T2indx)) ', ' noisetypes{nindx}])
        hold on
    end
end
xlim([0 1])
ylim([-0.05, 0.2])
xlabel('f_{IC}')
ylabel('Bias')
title(['\sigma_0 = ' num2str(sigma0)])
legend
saveas(fig, ['C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Figures\bias (Rice sigma0 ' num2str(sigma0rice) ') sigma0 ' num2str(sigma0) '.fig'])

% Variance figure
fig = figure;
for nindx = 1:length(noisetypes)
    for T2indx = 1:length(T2s)
        plot(fICs, squeeze(VarianceResults(nindx,T2indx,:)), noisemarkers{nindx}, color = T2colors{T2indx}, DisplayName=['T2: ' num2str(T2s(T2indx)) ', ' noisetypes{nindx}])
        hold on
    end
end
xlim([0 1])
ylim([0 0.03])
xlabel('f_{IC}')
ylabel('Variance')
title(['\sigma_0 = ' num2str(sigma0)])
legend
saveas(fig, ['C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\Figures\variance (Rice sigma0 ' num2str(sigma0rice) ') sigma0 ' num2str(sigma0) '.fig'])
%% Present as table