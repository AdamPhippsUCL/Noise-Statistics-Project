function [counts, bins, data] = simulateHistogram(opts)


arguments
    
    % OPTIONS
    
    % Signal
    opts.T2 = 100  % NaN 
    opts.fd 
    opts.TE 
    opts.Nav_ratio = 1

    % samples
    opts.nsample = 150

    % Noise
    opts.sigma0 = 0.05

    % Distribution type
    opts.disttype = 'Ratio'
    opts.zmin = 0
    opts.zmax = 2
    opts.dz = 0.005

    % Histogram bins
    opts.binmin = 0
    opts.binmax = 2
    opts.nbin = 100

end


%% Generate distribution

% b0 signal
b0signal = exp(-opts.TE/opts.T2);

% b signal
bsignal = opts.fd*b0signal;

switch opts.disttype

    case 'Ratio'

        [dist, signals] = RatioDistRician( ...
            b0signal,...
            bsignal,...
            opts.sigma0,...
            Nav_ratio = opts.Nav_ratio,...
            zmin = opts.zmin,...
            zmax = opts.zmax,...
            dz = opts.dz...
            );

    case 'Rice'

        disp('Need to write Rice distribution function')

end


%% Sample distribution

samples = zeros(opts.nsample, 1);

for sindx = 1:opts.nsample
    sample = sampleDistribution(dist, signals);
    samples(sindx) = sample;
end


%% Bin data into histogram

binedges = linspace(opts.binmin, opts.binmax, opts.nbin+1);

H = histogram(samples, binedges);

counts = H.Values;
bins = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;

data.binspacing = (opts.binmax-opts.binmin)/opts.nbin;


end