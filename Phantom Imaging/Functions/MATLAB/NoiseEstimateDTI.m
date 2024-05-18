% MATLAB function to use DTI images to estimate noise parameters (patch
% based approach)

function [sigma0, T2, resnorm] = NoiseEstimateDTI(Echo1fname, Echo2fname, opts)

arguments
    Echo1fname % DICOM fname of first echo image (TE1)
    Echo2fname % DICOM fname of first echo image (TE1)

    % OPTIONS
    opts.DATE = 'NA'

    % image info
    opts.multiframe = true
    opts.NSA = 4; % for corresponding diffusion images

    % Python script for saving images
    opts.pyscript = "C:\Users\adam\OneDrive - University College London\UCL PhD\Image-Processing\DICOM\Python\normaliseDiffusionImage.py";

    % Fitting
    opts.patchsize = 3
    opts.disttype = 'DTIRatio'
    opts.nbin = 20;
    opts.mask = [];

    % Folder to save mat images in
    opts.OutputFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Outputs\Noise Statistics"

end

%% STEPS

% 1. Read DICOM information for each echo image (Need their TE values)
% 2. Save/Load DICOM images as MAT files (CAREFUL WITH IMAGE ORIENTATION!)
% 2. Find ratio of echo images
% 3. Define mask/ROI for estimation
% 4. For each voxel, estimate sigma0 and T2 over surrounding patch (starting values from quick
% median calculation)
% 5. Return maps of sigma0 and T2


opts.ImageOutputFolder = [char(opts.OutputFolder) '/' opts.DATE '/Images'];

%% 1. Read DICOM information for each echo image (Need their TE values)

dinfo1 = dfparse(Echo1fname);
dinfo2 = dfparse(Echo2fname);

try
    TE1 = dinfo1(1).EffectiveEchoTime;
catch
    TE1 = dinfo1(1).EchoTime;
end

try
    TE2 = dinfo2(1).EffectiveEchoTime;
catch
    TE2 = dinfo2(1).EchoTime;
end

function [b, a] = swap(a, b)
end

% Make sure img2 has longer echo time
if TE2<TE1
    swap(dinfo1, dinfo2)
    swap(Echo1fname, Echo2fname)
    swap(TE1, TE2)
end



SeriesDescription1 = dinfo1(1).ProtocolName;
Ndirec1 = max([dinfo1.DiffGradOrientIdentifier]);
bval1 = dinfo1(1).DiffusionBValue;


SeriesDescription2 = dinfo2(1).ProtocolName;
Ndirec2 = max([dinfo2.DiffGradOrientIdentifier]);
bval2 = dinfo2(1).DiffusionBValue;
try
    TE2 = dinfo2(1).EffectiveEchoTime;
catch
    TE2 = dinfo2(1).EchoTime;
end


if Ndirec1 ~= Ndirec2
    error("Images have different number of DTI directions!");
else
    Ndirec = Ndirec1;
end



%% 2. Save/Load DICOM images as MAT files

Imagefnames = {Echo1fname, Echo2fname};
for indx = 1:length(Imagefnames)
    imagefname = Imagefnames{indx};
    pyrunfile( ...
        opts.pyscript, ...
        imagefname = imagefname, ...
        multiframe = opts.multiframe, ...
        returnimages = false,...
        saveimages = true,...
        datatype = 'mat',...
        outputfolder =  opts.ImageOutputFolder...
        );
end

%% OLD CODE WHICH WORKS
% % == Create stacked normalised image
% for direcIndx = 1:Ndirec
% 
% 
%     try
%         % Load first echo image
%         img1 = load([char(opts.ImageOutputFolder) '/mat/' char(SeriesDescription1) '/bval' num2str(bval1) '.0_DiffDirec' num2str(direcIndx) '.0.mat']).img;
% 
%         % Load second echo image
%         img2 = load([char(opts.ImageOutputFolder) '/mat/' char(SeriesDescription2) '/bval' num2str(bval1) '.0_DiffDirec' num2str(direcIndx) '.0.mat']).img;
%     catch
% 
%         % Load first echo image
%         img1 = load([char(opts.ImageOutputFolder) '/mat/' char(SeriesDescription1) '/bval' num2str(bval1) '.0_DiffDirec' num2str(direcIndx) '.mat']).img;
% 
%         % Load second echo image
%         img2 = load([char(opts.ImageOutputFolder) '/mat/' char(SeriesDescription2) '/bval' num2str(bval1) '.0_DiffDirec' num2str(direcIndx) '.mat']).img;
% 
%     end
%     % Remove zeros
%     img1(img1==0) = eps;
%     img2(img2==0) = eps;
% 
%     % imshow(squeeze( img1(8,:,:)) , [])
%     % imshow(squeeze( img1(8,:,:)) , [])
% 
%     % Calculate ratio
%     ratioimg = img2./img1;
% 
% 
%     if direcIndx == 1
%         RatioImageStack = zeros([size(ratioimg) Ndirec]);
%     end
% 
%     RatioImageStack(:,:,:,direcIndx) = ratioimg;
% 
% 
% end

%% NEW CODE

% == Create stacked normalised image
for direcIndx1 = 1:Ndirec
        
    try
        % Load first echo image
        img1 = load([char(opts.ImageOutputFolder) '/mat/' char(SeriesDescription1) '/bval' num2str(bval1) '.0_DiffDirec' num2str(direcIndx1) '.0.mat']).img;
    catch
        % Load first echo image
        img1 = load([char(opts.ImageOutputFolder) '/mat/' char(SeriesDescription1) '/bval' num2str(bval1) '.0_DiffDirec' num2str(direcIndx1) '.mat']).img;
    end
    
    for direcIndx2 = 1:Ndirec

        try
            % Load second echo image
            img2 = load([char(opts.ImageOutputFolder) '/mat/' char(SeriesDescription2) '/bval' num2str(bval1) '.0_DiffDirec' num2str(direcIndx2) '.0.mat']).img;
        catch
            % Load second echo image
            img2 = load([char(opts.ImageOutputFolder) '/mat/' char(SeriesDescription2) '/bval' num2str(bval1) '.0_DiffDirec' num2str(direcIndx2) '.mat']).img;
        end
    
        % Remove zeros
        img1(img1==0) = eps;
        img2(img2==0) = eps;
    
        % imshow(squeeze( img1(8,:,:)) , [])
        % imshow(squeeze( img1(8,:,:)) , [])
    
        % Calculate ratio
        ratioimg = img2./img1;
    
    
        % Configure array for ratio image stack
        if and(direcIndx1 == 1, direcIndx2 == 1)
            RatioImageStack = zeros([size(ratioimg) Ndirec]);
        end
    
        RatioImageStack(:,:,:, (direcIndx1-1)*Ndirec+direcIndx2 ) = ratioimg;

    end
end


%%
% Define blank arrays for T2 and sigma0 estimates and error
T2 = squeeze( zeros(size(RatioImageStack(:,:,:,1))) );
sigma0 = squeeze( zeros(size(RatioImageStack(:,:,:,1))) );
error = squeeze( zeros(size(RatioImageStack(:,:,:,1))) );

%% Define mask (for testing)

% Define slice and pixel indices in small ROI in capsule
zs = [1:size(ratioimg, 1)];
ys = [opts.patchsize:size(ratioimg, 2)-opts.patchsize];
xs = [opts.patchsize:size(ratioimg, 3)-opts.patchsize];

% Patch width
pw = ceil((opts.patchsize-1)/2);



%% Define histogram

for zindx = zs
    zindx
    for yindx = ys
        for xindx = xs
            
            patchvals = RatioImageStack(zindx, yindx-pw:yindx+pw, xindx-pw:xindx+pw, :);
            patchvals  = patchvals(:);
            
            % N Define outlier range
            iqr = (prctile(patchvals, 75) - prctile(patchvals, 25));
            up = median(patchvals) + 1.5*iqr;
            down = median(patchvals) - 1.5*iqr;
            
            % Remove outliers
            outlierbools = or(patchvals>up, patchvals<down);
            patchvals(outlierbools) = [];

            % == INITIAL GUESSES

            T2guess = -(TE2-TE1)/log(mean(patchvals));
            % Deal with unrealistic values
            T2guess(T2guess<25) = 25;
            T2guess(T2guess>250) = 250;


            rTE = (2*(TE1^2))/(TE1+TE2);

            sigma0guess = std(patchvals)/(sqrt(2)*exp(rTE/T2guess));
            sigma0guess(sigma0guess > 0.25) = 0.25;            
            sigma0guess(sigma0guess < 0.001) = 0.001;
           
            % % == USING RATIO DISTRIBUTION (SLOW)

            % % Define bin edges and centers
            % binmin = 0;
            % binmax = 2;
            % nbin = opts.nbin;
            % binedges = linspace(binmin, binmax, nbin+1);
            % bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
            % binspacing = bincenters(2)-bincenters(1);
            % 
            % 
            % counts = histcounts(patchvals, binedges);

            % [coeffs, resnorm] = fitDistToHist( ...
            %     counts, ...
            %     bincenters, ...
            %     fd = 1, ...
            %     TE = [TE1,TE2], ...
            %     disttype = opts.disttype, ...
            %     beta0guess = [sigma0guess, T2guess,1], ...
            %     boundstype = 'fixed',...
            %     lb = [0.001, 30, 0.1],...
            %     ub = [0.25, 400, 1]...
            %     );
            % 
            % sigma0fit = coeffs(1)
            % T2fit = coeffs(2)
            

            % == JUST USING INITIAL GUESSES (fast)
 
            sigma0fit = sigma0guess;
            T2fit = T2guess;
            resnorm = 1;



            
            % Append to arrays
            T2(zindx, yindx, xindx) = T2fit;
            sigma0(zindx, yindx, xindx) = sigma0fit;
            error(zindx, yindx, xindx) = resnorm;
            
            % % 
            % % Define histogram
            % % f = figure('visible','off');
            % f= figure;
            % H = histogram(patchvals, binedges);
            % hold on;
            % counts = H.Values;
            % 
            % % Generate fitted distribution
            % [dist, signals] = RatioDistRician(exp(-TE1/T2fit), exp(-TE2/T2fit), sigma0fit, Nav_ratio=1);
            % 
            % % Add to histogram
            % hold on
            % plot(signals, dist*sum(counts)*binspacing, '-*')
            % % xlim([0.5,1.5])
            % % T2fit
            % % sigma0fit
            % close(f);

        end
    end
end


%% Account for NSA
sigma0 = sigma0/sqrt(opts.NSA);

%% Smoothing

smoothsigma = 0.5;
% % Filter slices
for slice = zs
    T2(slice,:,:) = imgaussfilt(squeeze( T2(slice,:,:)), smoothsigma);
    sigma0(slice,:,:) = imgaussfilt(squeeze( sigma0(slice,:,:)), smoothsigma);
end

%% Save maps

% Save Output images as mat files
MapOutputFolder = [char(opts.OutputFolder) '/' opts.DATE '/Maps'];

save([MapOutputFolder '/T2.mat'], "T2")
save([MapOutputFolder '/sigma0.mat'], "sigma0")
save([MapOutputFolder '/img1.mat'], "img1")
save([MapOutputFolder '/img2.mat'], "img2")

% % % Filter slices
% slice = 5;
% % T2slice = squeeze(T2(slice,:,:));
% % sigma0slice = squeeze(sigma0(slice, :,:));
% T2slice = imgaussfilt(squeeze( T2(slice,:,:)), 1);
% sigma0slice = imgaussfilt(squeeze( sigma0(slice,:,:)), 1);
% 
% figure;
% imshow(squeeze(img1(slice,:,:)), [])
% colorbar
% figure;
% imshow(T2slice, [0, 200])
% colorbar
% figure;
% imshow(sigma0slice*0.5, [0.0 0.075])
% colorbar

end