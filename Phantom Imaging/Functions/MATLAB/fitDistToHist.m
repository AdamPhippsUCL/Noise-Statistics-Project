function [coeffs, resnorm] = fitDistToHist(counts, bincenters, params, opts)

arguments

    counts % Array of counts
    bincenters % Array of bin centers

    % Distribution parameters 
    params.sigma0
    params.T2
    params.fd
    params.TE = 50 % Specify at [TE1, TE2] for DTIRatio
    params.Nav_ratio = 1

    opts.disttype = 'DTIRatio'
    opts.beta0guess = [0.05, 250, 0.8]

    opts.boundstype = 'fixed';
    opts.lb = [0.001, 200, 0.001]
    opts.ub = [0.025, 1500, 1]
    opts.boundsfrac = 0.25

end


%% Steps

% 1. Normalise histogram counts
% 2. Define unknown parameter vector for fitting
% 3. Apply fitting

%% 1. Normalise histogram counts
counts = counts/sum(counts);

%% 2. Define unknown parameters for fitting

Nparam = 3;
if isfield(params, 'fd')
    if isnan(params.fd)
        params = rmfield(params, 'fd');
    else
        Nparam = Nparam-1; % Just fitting for sigma0 and T2; x = [sigma0, T2]
    end
end

if isfield(params, 'T2')
    if isnan(params.T2)
        params = rmfield(params, 'T2');
    else
        Nparam = Nparam-1; % Just fitting for sigma0; x = sigma0
    end
end



%% Fit

% Use nlm function: fitnlm

% Assign variable to base (for access in functions)
assignin('base', 'current_params', params)


% Define bounds
switch opts.boundstype

    case 'fixed'
        % Bounds set by user
        disp('');

    case 'dependent'
        opts.lb = (1-opts.boundsfrac)*opts.beta0guess;
        opts.ub = (1-opts.boundsfrac)*opts.beta0guess;


end



switch opts.disttype

    case 'Ratio'
        if Nparam == 1
            [coeffs, resnorm] = lsqcurvefit(@RatioDist1, opts.beta0guess(1), bincenters, counts, opts.lb(1), opts.ub(1) );   
        elseif Nparam == 2
            [coeffs, resnorm] = lsqcurvefit(@RatioDist2, opts.beta0guess(1:2), bincenters, counts, opts.lb(1:2), opts.ub(1:2) );
        elseif Nparam == 3
            [coeffs, resnorm] = lsqcurvefit(@RatioDist3, opts.beta0guess(1:3), bincenters, counts, opts.lb(1:3), opts.ub(1:3) );
        end

    case 'Rice'

        if Nparam == 1
            [coeffs, resnorm] = lsqcurvefit(@RiceDist1, opts.beta0guess(1), bincenters, counts, opts.lb(1), opts.ub(1:3) );  
        elseif Nparam == 2
            [coeffs, resnorm] = lsqcurvefit(@RiceDist2, opts.beta0guess(1:2), bincenters, counts, opts.lb(1:2), opts.ub(1:2) );   
        elseif Nparam == 3
            [coeffs, resnorm] = lsqcurvefit(@RiceDist3, opts.beta0guess(1:3), bincenters, counts, opts.lb(1:3), opts.ub(1:3) ); 
        end

    case 'DTIRatio'

        if Nparam == 1
            error(['Should be two parameter fitting, not ' Nparam ' parameters'])
        elseif Nparam == 2
            [coeffs, resnorm] = lsqcurvefit(@DTIRatioDist2, opts.beta0guess(1:2), bincenters, counts, opts.lb(1:2), opts.ub(1:2) );   
        elseif Nparam == 3
            error(['Should be two parameter fitting, not ' Nparam ' parameters'])
        end

    case 'DTIRice'

        if Nparam == 1
            error(['Should be two parameter fitting, not ' Nparam ' parameters'])
        elseif Nparam == 2
            [coeffs, resnorm] = lsqcurvefit(@DTIRiceDist2, opts.beta0guess(1:2), bincenters, counts, opts.lb(1:2), opts.ub(1:2) );   
        elseif Nparam == 3
            error(['Should be two parameter fitting, not ' Nparam ' parameters'])
        end


end

if Nparam == 1
    coeffs = [coeffs, params.T2, params.fd, params.TE];
elseif Nparam == 2
    coeffs = [coeffs, params.fd, params.TE];
elseif Nparam == 3
    coeffs = [coeffs, params.TE];
end


end


% Function to fit Ratio distribution

function binfreqs = RatioDist1(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % sigma0 value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
T2 = params.T2;
TE = params.TE;
fd = params.fd;
Nav_ratio = params.Nav_ratio;
sigma0 = b;

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RatioDistRician( ...
    b0signal, ...
    bsignal, ...
    sigma0, ...
    Nav_ratio=Nav_ratio, ...
    zmin = binmin, ...
    zmax = binmax, ...
    dz = binspacing,...
    ymin = binmin,...
    ymax = binmax,...
    dy = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end

function binfreqs = RatioDist2(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2] value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
TE = params.TE;
fd = params.fd;
Nav_ratio = params.Nav_ratio;
sigma0 = b(1);
T2 = b(2);

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RatioDistRician( ...
    b0signal, ...
    bsignal, ...
    sigma0, ...
    Nav_ratio=Nav_ratio, ...
    zmin = binmin, ...
    zmax = binmax, ...
    dz = binspacing,...
    ymin = binmin,...
    ymax = binmax,...
    dy = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));


end

function binfreqs = RatioDist3(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2] value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
TE = params.TE;
Nav_ratio = params.Nav_ratio;
sigma0 = b(1);
T2 = b(2);
fd = b(3);

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RatioDistRician( ...
    b0signal, ...
    bsignal, ...
    sigma0, ...
    Nav_ratio=Nav_ratio, ...
    zmin = binmin, ...
    zmax = binmax, ...
    dz = binspacing,...
    ymin = binmin,...
    ymax = binmax,...
    dy = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));


end


% Function to fit Rice distribution

function binfreqs = RiceDist1(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % sigma0 value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
T2 = params.T2;
TE = params.TE;
fd = params.fd;
Nav_ratio = params.Nav_ratio;
sigma0 = b;

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RiceDist(b0signal, bsignal, sigma0, Nav_ratio=Nav_ratio, zmin = binmin, zmax = binmax, dz = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end

function binfreqs = RiceDist2(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2] values
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
TE = params.TE;
fd = params.fd;
Nav_ratio = params.Nav_ratio;
sigma0 = b(1);
T2 = b(2);

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RiceDist(b0signal, bsignal, sigma0, Nav_ratio=Nav_ratio, zmin = binmin, zmax = binmax, dz = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end

function binfreqs = RiceDist3(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2, fd] values
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
TE = params.TE;
Nav_ratio = params.Nav_ratio;
sigma0 = b(1);
T2 = b(2);
fd = b(3);

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RiceDist(b0signal, bsignal, sigma0, Nav_ratio=Nav_ratio, zmin = binmin, zmax = binmax, dz = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end


% Function to fit DTI Ratio distribution

function binfreqs = DTIRatioDist2(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2] values
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
TE1 = params.TE(1);
TE2 = params.TE(2);
fd = params.fd;
Nav_ratio = params.Nav_ratio;
sigma0 = b(1);
T2 = b(2);

% Echo 1 signal
Echo1signal = exp(-TE1/T2);

% Echo 2 signal
Echo2signal = exp(-TE2/T2);


% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RatioDistRician(Echo1signal, Echo2signal, sigma0, Nav_ratio=Nav_ratio, zmin = binmin, zmax = binmax, dz = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end


% Function to fit DTI Rice distribution

function binfreqs = DTIRiceDist2(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2] values
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
TE1 = params.TE(1);
TE2 = params.TE(2);
fd = params.fd;
Nav_ratio = params.Nav_ratio;
sigma0 = b(1);
T2 = b(2);

% Echo 1 signal
Echo1signal = exp(-TE1/T2);

% Echo 2 signal
Echo2signal = exp(-TE2/T2);


% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RiceDist(Echo1signal, Echo2signal, sigma0, Nav_ratio=Nav_ratio, zmin = binmin, zmax = binmax, dz = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end

