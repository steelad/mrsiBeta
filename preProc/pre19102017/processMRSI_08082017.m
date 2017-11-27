function [dw_mean_hlsvd,img_water,meta_mask,mvmt,img1] = ...
    processMRSI_08082017(subj,GF,LB,t1Thresh,mot,PATH,outSuffix)
%% preprocessing script standardized through August 8 2017
% inputs:
%  - subj:          subject number
%  - GF:            gaussian filter applied to timeseries
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
%subjs = {'3','4','5','6'};

%for INDEX = 1:5

%    clearvars -except INDEX subjs

%subj = 'HV_101'
%disp(strcat('subj00',subjs{INDEX}));
%subj = strcat('subj00',subjs{INDEX});
if isempty(PATH) == 0
    cd(PATH);
end
outName = strcat(subj,'_thruHSVD_limVOI_',outSuffix,'.mat');
%%
addpath('/home/fs0/asteel//scratch/uzayScripts/scripts_templates/scripts/mrsiBeta')
addpath('/home/fs0/asteel//scratch/uzayScripts/irt')
%addpath('/home/fs0/asteel//scratch/uzayScripts/subject_1')
addpath(genpath('/home/fs0/asteel//scratch/uzayScripts/SPID/'));
%%
% try
%     cc1=load(strcat(subj,'_01_metab_water_dw_limitedvoi.mat'));
%     cc2=load(strcat(subj,'_02_metab_water_dw_limitedvoi.mat'));
%     cc3=load(strcat(subj,'_03_metab_water_dw_limitedvoi.mat'));
% catch
% 
%     cc1=load(strcat(subj,'_01_metab_water_dw_limtedvoi.mat'));
%     cc2=load(strcat(subj,'_02_metab_water_dw_limtedvoi.mat'));
%     cc3=load(strcat(subj,'_03_metab_water_dw_limitedvoi.mat'));
% end
cc= struct;
imgFiles = dir([subj '_metab_water*']);
imgs = cell(length(imgFiles),1);
%%
for i = 1:length(imgFiles)
    imgs{i} = load(imgFiles(i).name);
end

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
%%
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
%% choose mask limits
if isempty(t1Thresh)
    accept = 0;
    redraw = 0;
    figure
    colormap gray
    imagesc(max(temp,[],3))
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
temp2 = max(temp,[],3);
meta_mask = zeros(size(temp2));
meta_mask(temp2 > lowerLim) = 1;
meta_mask2 = meta_mask;
[rows,cols,~] = find(meta_mask);
meta_mask2(min(rows)-2:max(rows)+2,min(cols)-2:max(cols)+2) = 1; 
meta_mask2(min(rows)-2:max(rows)+2,min(cols)-2:max(cols)+2) = 1; 
meta_mask2(min(rows):min(rows)-2,min(cols)-2:max(cols):max(cols)+2) = 1; 
meta_mask2(max(rows):max(rows)+2,min(cols)-2:max(cols):max(cols)+2) = 1; 

%%
if isempty(t1Thresh)
%% draw and confirm
    close all
    figure
    subplot(1,2,1); imagesc(meta_mask2); 
    axis square
    subplot(1,2,2); imagesc(max(temp,[],3))
    input('++ Press enter to continue with this mask')
    axis square
    clear accept currLim
    close all
end

%%
for i = 1:size(temp,3)
    temp2 = temp(:,:,i);
    temp2(meta_mask2 == 0) = 0;
    temp(:,:,i) = temp2;
end
%%
if mot == 1
    disp('++ Start motion correction with rigid body registration')
    %% register images
        transMat = cell(1,numImgs);
        [optim,met] = imregconfig('monomodal');
        prop = imref2d(size(temp(:,:,1)));
        water_img1 = zeros(size(temp));
        mvmt = zeros(3,3,numImgs);
    for i = 1:numImgs
        displacementOutmatrix = imregtform(temp(:,:,i),temp(:,:,1),...
            'rigid',optim,met);
        transMat{i} = displacementOutmatrix;
        [water_img1(:,:,i),~] = imwarp(temp(:,:,i),...
            prop,displacementOutmatrix,'outputview',prop,'interp','nearest');
        mvmt(:,:,i) = transMat{i}.T;
    end
%%
mvmt = reshape(mvmt,[],10);
figure
subplot(2,1,1); title(['rotation (q)'])
plot(sin(mvmt(2,:)),'b')
subplot(2,1,2); title(['translation (total, vxl)'])
plot(mvmt(3,:).^2.*mvmt(6,:).^2,'r')
saveas(gcf,'motion.png','png')
close gcf
temp = img_dw_90_water;
for i = 1:numImgs
    disp(['++ Applying registration: ' num2str(i)])
    for j = 1:numTPs
        img_dw_90(:,:,j,i) = imwarp(img_dw_90(:,:,j,i),...
        prop,transMat{i},'outputview',prop,'interp','nearest');
        img_dw_90_water(:,:,j,i) = imwarp(img_dw_90_water(:,:,j,i),...
        prop,transMat{i},'outputview',prop,'interp','nearest');
    end
end
makeImgTimelapse(temp,img_dw_90_water,'mvmt.mp4')
close all
else
    mvmt = [];
%%
end

%%
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
dw=1/1250;
tt=0:dw:dw*(1024-1);
%GF=0.2;
%LB=0;
LB1=LB;
disp('++ Begin frequency cross correlation.')
disp('++ This may take a while')
for kk=1:48
for jj=1:48
for zz=1:numImgs*2
if (mod(zz,2))
    cc1=[conj(squeeze(abc_dw(kk,jj,:,zz))); zeros(512,1)];
         aa(:,zz)=cc1.*(exp(-tt*pi*LB-tt.^2/(GF^2)))';
else
    cc1=[conj(squeeze(abc_dw(kk,jj,:,zz))); zeros(512,1)];
    
  aa(:,zz)=cc1.*(exp(-tt*pi*LB1-tt.^2/(GF^2)))';
end
 end
 bb = freq_xcorr_uzay(aa ,1250,123.2);%works better

          bb=phase_max_uzay(bb);
  img1(kk,jj,:,:) =  bb;      
end

end
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
saveas(gcf,'subtraction_img1.png','png')
close gcf
end;
%%
% c = 0;
% %%
% for i = 1:2:6
%     c = c+1;
    dw_mean(:,:,:)=mean(-squeeze(img1(:,:,:,1:2:end)),4)+...
        mean(1.005*squeeze(img1(:,:,:,2:2:end)),4);
% end
%%
figure
plot((real((fft(squeeze(conj(dw_mean(19,22,:))))))));
saveas(gcf,'example_19_22.png','png')
close gcf
%%
c = 1
ppmrange = [0.0 4.2];
 %ppmrange = [6 9.9];
%ppmrange = [0.5 4.2];%5_11_bg
step = 1/1250*1000;
frequency= 123.2*1000;
fsampl=1250
ndp=1024
xmin = (ppmrange(1)-4.7)*frequency/10^6;%in kHz
xmax = (ppmrange(2)-4.7)*frequency/10^6;%in kHz
t = [0:step:(ndp-1)*step];
freqrange1 = [-(fsampl/2-1),xmin*1000]; %in Hz
freqrange2 = [xmax*1000,(fsampl/2-1)]; %in Hz
freqrange = [freqrange1;freqrange2];
MM = 1024/2;
M=30;
begin=0;
boundL = 0; %in ppm
boundH = 9; % in ppm
%load('subject_6_9_11_vapor_bg.spa','-mat')
dw=1/1250;
tt=0:dw:dw*(1024-1);
%GF=10000;
%LB=0;
disp('++ Begin HLSVD. This will take a while...')
%for c = 1:3
for i=1:48
    for k=1:48
[dw_mean_hlsvd(i,k,:,c),~,~,~,~,~] = HLSVDPROrec((squeeze(dw_mean(i,k,:,c))),M,fsampl,t/1000,MM,freqrange);
    end
end
%end
%%
figure
plot((real((fft((squeeze(conj(dw_mean_hlsvd(22,19,:)))))))));
saveas(gcf,'post_hlsvd_19_22.png','png')
close gcf
%%
img_water = mean(img_dw_90_water,4);
%%
save(outName,'dw_mean_hlsvd','img_water','meta_mask','mvmt','img1','-mat');
%end
if mot == 1
    genspa(subj,outName,['LB' num2str(LB) '_wMot'])
else
    genspa(subj,outName,['LB' num2str(LB) '_noMot'])
end
