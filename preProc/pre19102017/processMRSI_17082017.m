function [dw_mean_hlsvd,img_water,meta_mask,mvmt,img1,outName] = ...
    processMRSI_17082017(subj,GF,LB,ave,t1Thresh,mot,PATH,outSuffix)
%% preprocessing script standardized through August 8 2017
% inputs:
%  - subj:          subject number
%  - GF:            gaussian filter applied to timeseries
%  - avg:           Do some averaging
%  - LB:            line broadening level
%  - t1Thresh:      Draw mask based on T1 thresholded value
%  - mot:           Apply motion correction and output motion value
%  - PATH:          Path to data
%  - outSuffix:     suffix for the thruHSVD file
% outputs:
%  - dw_mean_hlsvd: Fully cleaned timeseries for spectra analysis
%  - img_water:     Water image
%  - meta_mask:     Metabolite mask drawn via T1 thresholding
%  - mvmt:          Movement for each time point
%  - img1:          timeseries before smoothing & hlsvd, etc.
%%

if isempty(PATH) == 0
    cd(PATH);
end
outName = strcat(subj,'_thruHSVD_',outSuffix,'.mat');
%%
addpath('/home/fs0/asteel//scratch/uzayScripts/scripts_templates/scripts/mrsiBeta')
addpath('/home/fs0/asteel//scratch/uzayScripts/irt')
addpath(genpath('/home/fs0/asteel//scratch/uzayScripts/SPID/'));
cc= struct;
imgFiles = dir([subj '*_metab_water*']);
imgs = cell(length(imgFiles),1);
%%
for i = 1:length(imgFiles)
    disp(['++ Loading image ' num2str(i)])
    imgs{i} = load(imgFiles(i).name);
end
%%
%%
%%
for i = 1:length(imgs)
    if i == 1
        cc.img_metab = imgs{1}.img_metab;
        cc.img_water = imgs{1}.img_water;
    else
        cc.img_metab = cat(4,cc.img_metab,imgs{i}.img_metab);
        cc.img_water = cat(4,cc.img_water,imgs{i}.img_water);
    end
end
%%

%%
numImgs = length(imgs);
numTPs = size(cc.img_metab,3);
temp = cc.img_metab(:,:,1,1);
%%
img_dw_90 = zeros([size(fliplr(imrotate(squeeze(temp),-90))) numTPs numImgs*2]);
img_dw_90_water = zeros([size(fliplr(imrotate(squeeze(temp),-90))) numTPs numImgs]);
disp('++ Doing some flipping and rotating for display purposes')
for kk=1:numImgs*2
    for jj=1:numTPs

        img_dw_90(:,:,jj,kk)=fliplr(imrotate(squeeze(cc.img_metab(:,:,jj,kk)),-90));

    end
end

for kk=1:numImgs
    for jj=1:numTPs

        img_dw_90_water(:,:,jj,kk)=fliplr(imrotate(squeeze(cc.img_water(:,:,jj,kk)),-90));

    end
end
%%
temp = squeeze(sum(abs(img_dw_90_water),3));
disp('++ Creating metab masks')
%% choose mask limits
if ~ischar(t1Thresh)
    [meta_mask] = threshMetaMask(temp,t1Thresh);
    lipid_mask = [];
else
    disp('++ Whole brain requested. Draw masks by hand later.')
end

%%

%%
if mot ~= 0
    disp('++ Doing motion correction')
    disp('++ This might take a while')
    [img_dw_90_water,img_dw_90,mvmt] = motMRSI(temp,img_dw_90_water,img_dw_90,mot);
else
    disp('Skipping motion correction')
    mvmt = [];
end
%%
abc_dw = zeros(size(img_dw_90));
%%
disp('++ Eddy correct')
for kk=1:48
    for ii=1:48
        for ll=1:numImgs*2
            if mod(ll,2) == 1
                wInd = (ll+1)/2;
                abc_dw(kk,ii,:,ll)=spec_ecc_cor(conj(squeeze(img_dw_90(kk,ii,:,ll))),...
                    conj(squeeze(img_dw_90_water(kk,ii,:,wInd))),0,0);
            else
                wInd = ll/2;
                abc_dw(kk,ii,:,ll)=spec_ecc_cor(conj(squeeze(img_dw_90(kk,ii,:,ll))),...
                    conj(squeeze(img_dw_90_water(kk,ii,:,wInd))),0,0);
            end
        end
    end
end


%%

disp('++ Begin frequency cross correlation.')
disp('++')
disp('++ This may take a while')
%%
[img1] = do_phaseXcorr(abc_dw,LB,GF);
%%
if mod(size(img1,4),4) == 0
    c = 0;
        figure('position',[0 0 300 2000])
    for i = 2:2:size(img1,4)
        c = c+1;
        subplot(size(img1,4)/4,2,c)
        ind = i;
        plot((real((fft(squeeze(conj(img1(24,24,:,ind)))...
            -1.005*squeeze(conj(img1(24,24,:,ind-1))))))));
        ylim([-10 300])
        xlim([50 300])
    end
    
    saveas(gcf,['subtraction_img1_' outName '.png'],'png')
    close gcf
end;
%%
% c = 0;
% %%
% for i = 1:2:6
%     c = c+1;
if ave == 1
    dw_mean(:,:,:)=mean(-squeeze(img1(:,:,:,1:2:end)),4)+...
        mean(1.005*squeeze(img1(:,:,:,2:2:end)),4);
    figure
    plot((real((fft(squeeze(conj(dw_mean(19,22,:,1))))))));
    saveas(gcf,['example_19_22_' outName '.png'],'png')
    close gcf

else
    dw_mean(:,:,:,:)=-squeeze(img1(:,:,:,1:2:end))+...
        1.005*squeeze(img1(:,:,:,2:2:end));
end
% end
%%

%%

ppmrange = [0.0 4.2];
step = 1/1250*1000;
frequency= 123.2*1000;
fsampl=1250;
ndp=1024;
xmin = (ppmrange(1)-4.7)*frequency/10^6;%in kHz
xmax = (ppmrange(2)-4.7)*frequency/10^6;%in kHz
t = 0:step:(ndp-1)*step;
freqrange1 = [-(fsampl/2-1),xmin*1000]; %in Hz
freqrange2 = [xmax*1000,(fsampl/2-1)]; %in Hz
freqrange = [freqrange1;freqrange2];
MM = 1024/2;
M=30;
disp('++ Begin HLSVD. This will take a while...')
dw_mean_hlsvd = zeros(size(dw_mean));
for c = 1:numImgs
    for i=1:48
        for k=1:48
            [dw_mean_hlsvd(i,k,:,c),~,~,~,~,~] = HLSVDPROrec((squeeze(dw_mean(i,k,:,c))),M,fsampl,t/1000,MM,freqrange);
        end
    end
    disp(['++ Finished HLSVD: ' num2str(c) ' of ' num2str(numImgs')])
end
%end
%%
if ave == 1
    figure
    plot((real((fft((squeeze(conj(dw_mean_hlsvd(22,19,:)))))))));
    saveas(gcf,['post_hlsvd_19_22_' outName '.png'],'png')
    close gcf
end
%%
if ave == 1
    img_water = mean(img_dw_90_water,4);
else
    img_water = img_dw_90_water;
end
%%
if ~ischar(t1Thresh)
    save(outName,'dw_mean_hlsvd','img_water','meta_mask','mvmt','img1','-mat');
else
    save(outName,'dw_mean_hlsvd','img_water','mvmt','img1','-mat');
end
    
%%
if ischar(t1Thresh)
    [dw_mean_hlsvd,preLipidHLSVD,meta_mask,~,~] =...
    doWholeBrainLipid(subj,dw_mean_hlsvd,img_water,outName)
    
    %%
end
%end
if ~ischar(t1Thresh)
    if mot == 1
        genspa(subj,outName,['LB' num2str(LB) '_GF_' num2str(GF) '_wMot_av' num2str(ave)])
    elseif mot == 2
        genspa(subj,outName,['LB' num2str(LB)  '_GF_' num2str(GF) '_wMotLin_av' num2str(ave)])
    else
        genspa(subj,outName,['LB' num2str(LB)  '_GF_' num2str(GF) '_noMot_av' num2str(ave)])
    end
else
    if ~strcmp(t1Thresh,'tum')
        if mot == 1
            genspa(subj,outName,['LB' num2str(LB) '_wMot'])
        elseif mot == 2
            genspa(subj,outName,['LB' num2str(LB) '_wMotLin'])
        else
            genspa(subj,outName,['LB' num2str(LB) '_noMot'])
        end
    else
        genspa_TE110(subj,outName,['LB' num2str(LB)])
    end
end
   disp('++ .Done. ++')     
