% MATLAB Function to apply MLP fitting
function [fIC, fEES, fVASC, R, rmse] = ns_verdict_MLP_fit(schemename, modeltype, Y, opts)

arguments
    schemename % Name of scheme (char array)
    modeltype % Name of model (char array)
    Y % signal data [..., nscheme]

    % == OPTIONS
   
    % Noise level used in model training (either single values or image
    % arrays)
    opts.noisetype = 'Rice'
    opts.sigma0train = 0.05
    opts.T2train = 100

    % T2 and sigma0 values of trained MLPs
    opts.possibleT2s = [100];
    opts.possiblesigma0s = [0.05]

    % mask
    opts.mask = []

    % Scheme structure for checking
    opts.scheme = []

    % == Folders
    opts.FolderPath = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations'
    % Schemes
    opts.schemesfolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\General Code\Model Fitting\MLP\My Tests\Schemes'

    % Python MLP folder
    opts.pythonfolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\General Code\Model Fitting\MLP\My Tests\Python'

    opts.modelsfolder = 'C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\VERDICT Simulations\MLP Models'
end

% Scheme folder
schemesfolder = opts.schemesfolder;%[opts.FolderPath '\Schemes'];

% Load scheme
load([schemesfolder '/' schemename '.mat']);


% Check scheme agrees
if ~isempty(opts.scheme)

    % Check loaded and expected schemes agree 
    if ~isequal([scheme.bval], [opts.scheme.bval]) && isequal([scheme.delta], [opts.scheme.delta])
        error("schemes don't agree!")
    end

end



% REMOVE NANS AND INFINTIES
Y(isnan(Y)) = 0;
Y(isinf(Y)) = 0;


% Scheme, Data, and Image sizes
nscheme = length(scheme) ;
szY = size(Y) ;
szmap = szY(1:end-1) ;


% Apply mask
if ~isempty(opts.mask)
    % Check mask size
    if szmap == size(opts.mask)    
        Y = Y.*opts.mask;
    else
        disp('Mask size mismatch')
    end
end

% Flatten
Y = reshape(Y,[prod(szmap) nscheme]) ;


%% Create masks for different T2 and sigma0 values

% If single values given as input
if isscalar(opts.T2train)
    opts.T2train = opts.T2train*ones(szmap);
end
if isscalar(opts.sigma0train)
    opts.sigma0train = opts.sigma0train*ones(szmap);
end

% Round values in T2train array to possible T2 value
roundedT2s = roundtowardvec(opts.T2train, opts.possibleT2s);

% Round values in T2train array to possible T2 value
roundedsigma0s = roundtowardvec(opts.sigma0train, opts.possiblesigma0s);


%% Run MLP fitting python script

pyfname = [opts.pythonfolder '/callMLP.py' ];

indx = 0;
for T2 = opts.possibleT2s
    for sigma0 = opts.possiblesigma0s
        indx = indx+1;

        thisschemename = [schemename '/' opts.noisetype '/T2_' num2str(T2) '/sigma_' num2str(sigma0)];
        % thisschemename = schemename;

        % Load some META data about MLP model
        load([opts.FolderPath '/MLP Models/' modeltype '/' thisschemename '/Meta.mat'])
        
        % Radii used in fitting
        Rs = FitMeta.TrainingDataMeta.Rs;

        % Create mask for signals
        mask = and(roundedT2s == T2, roundedsigma0s == sigma0);
        mask = mask(:);
        Ymask = repmat(mask, 1, nscheme);

        thisY = Y.*mask;
      
        thisx = pyrunfile( ...
            pyfname, ...
            'x', ...
            signals = thisY, ...
            modeltype = modeltype, ...
            schemename = schemename, ...
            noisetype = opts.noisetype,...
            sigma0train = sigma0, ...
            T2train = T2,...
            modelsfolder = opts.modelsfolder);

        thisx = double(thisx);

        if indx == 1
            Xmask = repmat(mask, 1, size(thisx,2));
            x = thisx.*Xmask;
        else
            Xmask = repmat(mask, 1, size(thisx,2));
            x = x + thisx.*Xmask;
        end

    end
end

%% Reformat results

switch modeltype

    case 'Original VERDICT'
    
        fIC = sum(x(:,1:end-2),2);
        fEES = x(:,end-1);
        fVASC = x(:,end);
        R = sum( Rs.*x(:,1:end-2), 2); 

    case 'No VASC VERDICT'
    
        fIC = sum(x(:,1:end-1),2);
        fEES = x(:,end);
        fVASC = zeros(size(fEES));
        R = sum( Rs.*x(:,1:end-1), 2); 

    case 'RDI'
    
        muRs = FitMeta.TrainingDataMeta.muRs;
        fIC = sum(x(:,1:end-1),2);
        fEES = x(:,end);
        fVASC = zeros(size(fEES));
        R = sum( muRs.*x(:,1:end-1), 2); 

    case 'RDI v1.3'
        fIC = x(:,1);
        R = x(:,2);
        fEES = x(:,3);
        fVASC = zeros(size(fEES));

    case 'RDI v1.4'

        fIC = x(:,1);
        R = x(:,2);
        fEES = x(:,4);
        fVASC = zeros(size(fEES));
end


%% Reshape results
if length(szmap)>1
    fIC = reshape(fIC,szmap) ;
    fEES = reshape(fEES,szmap) ;
    fVASC = reshape(fVASC,szmap);
    R = reshape(R, szmap);
end

rmse = zeros(size(fIC));
end





function X=roundtowardvec(X,roundvec,type)
%function newnums=roundtowardvec(X,[roundvec],[type])
%
% This function rounds number(s) toward given values. If more than one
% number is given to round, it will return the matrix with each rounded
% value, otherwise it will return the single rounded value. It will ignore
% NaNs and return them back with NaNs.
%
% Inputs: X: the number(s) that you want rounded
%
%         roundvec:(opt) the values to round X to. If none given, it will
%           default to -inf:1:inf (and use the built in functions).
%
%         type:(opt) specifies which kind of rounding you want
%           the function to use.
%
%           Choices are: 'round' - round to nearest value
%                        'floor' - round toward -Inf
%                        'ceil'  - round toward Inf
%                        'fix'   - round toward 0
%                        'away'  - round away from 0 (ceil if positive and floor if negative)
%                     (see help files for more clarity)
%
%           If no type is given, the function will default to rounding to
%           the nearest value.
%
% Outputs: newnums: rounded values, in same shape as X input matrix
%          indices: indices of rounded values in roundvec
if nargin==0
	help roundtowardvec; %if nothing given, tell what to give
	return
elseif isempty(X)
	%if given empty, return empty without going through whole script
	return
end
if nargout>1
	error('Too many output variables are given');
end
if ~exist('type','var') || isempty(type)
	type='round';  %%round to nearest value if not specified
end
if ~exist('roundvec','var') || isempty(roundvec) || all(isnan(roundvec))
	if strcmpi(type,'round')
		%to nearest integer
		X=round(X);
	elseif strcmpi(type,'away')
		%nearest integer away from 0
		X=ceil(abs(X)).*sign(X);
	elseif strcmpi(type,'fix')
		%nearest integer toward 0
		X=fix(X);
	elseif strcmpi(type,'floor')
		%nearest integer toward -inf
		X=floor(X);
	elseif strcmpi(type,'ceil')
		%nearest integer toward inf
		X=ceil(X);
	else
		error('%sRound type not recognized. Options are:\n''round'' - round to nearest value\n''floor'' - round toward -Inf\n''ceil''  - round toward Inf\n''fix''   - round toward 0\n''away''  - round away from 0','')
	end
else
	%Ignore nan in roundvec
	roundvec(isnan(roundvec))=[];
	
	%Record which values are nan to ignore
	Xnan=isnan(X);
	
	%Hold onto size for returning value
	sz=size(X);
	
	%Calculate differences
	X=X(:);
	roundvec=roundvec(:)';
	diffs=bsxfun(@minus,X,roundvec);
	
	if strcmpi(type,'round') %to nearest value
		[~,inds]=min(abs(diffs),[],2);
		X=roundvec(inds);
	elseif strcmpi(type,'fix') %to nearest value toward 0
		
		iless=X<0;
		X(iless)=roundtowardvec(X(iless),roundvec,'ceil');
		X(~iless)=roundtowardvec(X(~iless),roundvec,'floor');
	elseif strcmpi(type,'ceil') %nearest value toward inf
		diffs(diffs>0)=nan;
		[~,inds]=min(abs(diffs),[],2);
		
		i_inf=X>max(roundvec);
		X=roundvec(inds);
		X(i_inf)=inf;
	elseif strcmpi(type,'floor') %nearest value toward -inf
		diffs(diffs<0)=nan;
		[~,inds]=min(abs(diffs),[],2);
		
		i_inf=X<min(roundvec);
		X=roundvec(inds);
		X(i_inf)=-inf;
	elseif strcmpi(type,'away') %nearest value away from 0
		
		iless=X<0;
		X(~iless)=roundtowardvec(X(~iless),roundvec,'ceil');
		X(iless)=roundtowardvec(X(iless),roundvec,'floor');
	else
		error('%sRound type not recognized. Options are:\n''round'' - round to nearest value\n''floor'' - round toward -Inf\n''ceil''  - round toward Inf\n''fix''   - round toward 0\n''away''  - round away from 0','')
	end
	
	%Return to output side
	X=reshape(X(:),sz);
	
	%Ignore nan in input dataset
	X(Xnan)=nan;
end
end