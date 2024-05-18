n=size(FittingVariances,3);
nprot = size(FittingBiases,2);
maxbias = max(FittingBiases(:));
minbias = min(FittingBiases(:));
maxvar = max(FittingVariances(:));

for protindx = 1:nprot

    figure;
    scatter(ones(n,1), squeeze(FittingBiases(1,protindx,:)))
    hold on
    scatter(2*ones(n,1), squeeze(FittingBiases(2,protindx,:)))    
    boxplot(transpose(squeeze(FittingBiases(:,protindx,:)) ))
    ylim([1.1*minbias, 1.1*maxbias])
    title(['bias' num2str(protindx)])

    figure;
    scatter(ones(n,1), squeeze(FittingVariances(1,protindx,:)), DisplayName = )
    hold on
    scatter(2*ones(n,1), squeeze(FittingVariances(2,protindx,:)))    
    boxplot(transpose(squeeze(FittingVariances(:,protindx,:)) ))
    ylim([0, 1.1*maxvar])
    title(['variance' num2str(protindx)])
end