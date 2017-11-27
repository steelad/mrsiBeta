function [meta_mask,meta_mask2] = threshMetaMask(imgWaterSum,t1Thresh)
    if isempty(t1Thresh)
    accept = 0;
    figure
    colormap gray
    imagesc(max(imgWaterSum,[],3))
    axis square
    curLim = get(gca,'clim');
    %%
    while accept == 0
        currkey = input(['++ Enter new clim, or nothing if happy. \n' ...
        'Press enter to continue: ']);
        if isempty(currkey)
            lowerLim = min(get(gca,'clim'));
            accept = 1;
        else
            caxis([currkey max(curLim)])

            axis square
            drawnow
        end
    end
    else
        lowerLim = t1Thresh;
    end
%% make meta mask
imgWaterSum2 = max(imgWaterSum,[],3);
meta_mask = zeros(size(imgWaterSum2));
meta_mask(imgWaterSum2 > lowerLim) = 1;
meta_mask2 = meta_mask;
% [rows,cols,~] = find(meta_mask);
% meta_mask2(min(rows)-2:max(rows)+2,min(cols)-2:max(cols)+2) = 1; 
% meta_mask2(min(rows)-2:max(rows)+2,min(cols)-2:max(cols)+2) = 1; 
% meta_mask2(min(rows):min(rows)-2,min(cols)-2:max(cols):max(cols)+2) = 1; 
% meta_mask2(max(rows):max(rows)+2,min(cols)-2:max(cols):max(cols)+2) = 1;
meta_mask2 = ones(size(meta_mask));

%%
if isempty(t1Thresh)
%% draw and confirm
    close all
    figure
    subplot(1,2,1); imagesc(meta_mask2); 
    axis square
    subplot(1,2,2); imagesc(max(imgWaterSum,[],3))
    input('++ Press enter to continue with this mask')
    axis square
    clear accept currLim
    close all
end