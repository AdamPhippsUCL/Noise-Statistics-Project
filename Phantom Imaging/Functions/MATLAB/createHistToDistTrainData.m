function createHistToDistTrainData(opts)

%% NEED TO THINK ABOUT THIS BECAUSE fd, TE, and T2 can vary!!

arguments

    % Whether T2
    opts.T2fit = false

    opts.sigma0range = [0.005, 0.25]
    opts.T2range = [25, 400]

    opts.savedata = true
    opts.savefolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\VERDICT Screening\Code\VERDICT-Screening\Noise Statistics\MLP\Training Data\HistToDist"

end



end