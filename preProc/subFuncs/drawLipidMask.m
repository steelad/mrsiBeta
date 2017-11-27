function [BW, BW1,lipid_mask] = drawLipidMask(imgWater,imgMetab)

%%
mvon = 0;
while mvon ~= 1
    %%
    figure; title('++ Draw Brain','Fontsize',16)
    ax1 = subplot(1,2,1);
    %%
    imagesc(sum(abs(imgWater(:,:,:)),3));
    axis square
tTEMP = fftshift(fft(conj(squeeze(imgMetab(24,24,:,1)))));
temp = zeros(size(imgMetab,1),size(imgMetab,2),length(tTEMP),1);
        for i=1:48
        for j=1:48
            temp(i,j,:)=fftshift(fft(conj(squeeze(imgMetab(i,j,:,1))))); % FFT to get into frequency domain, 
                                                                            % necessary for lipid 
        end
        end
        subplot(1,2,2)
        imagesc(sum(abs(temp(:,:,:)),3));
        axis square

    %%
    h = imfreehand(ax1);
    %%
    BW = createMask(h); %% brain
    %%
    figure; 
    imagesc(BW); 
    axis square
    mvon = input('++ Press 1 to move on. Press enter to redraw.');
end
%%
%%
%%
close all
mvon = 0;
while mvon ~=1
%%
    figure; title('++ Draw Lipid Extent','fontsize',16)
    %%
    imagesc(sum(abs(imgWater(:,:,:)),3));
    axis square
    h = imfreehand();
    BW1 = createMask(h); %% outer layer
    %%
    lipid_mask=BW1-BW;  %% subtraction
    %%
    close all
    %%
    figure
    imagesc(lipid_mask)
    axis square
    mvon = input('++ Press enter again to move on with this lipid mask');
end
    close all
%%

%%
