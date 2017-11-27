function [dw_mean_hlsvd,img_water,meta_mask] = ...
    procWithExcl(subj,insuff,outSuffix,excl,PATH)
cd([PATH '/' subj])
pwd
load([subj insuff '.mat'])
%%
EXCL = excl.*2;
evExcl = 2:2:size(img1,4); oddExcl = 1:2:size(img1,4);
evExcl = setdiff(evExcl,EXCL); oddExcl = setdiff(oddExcl,EXCL-1);
%%
subtr = -squeeze(img1(:,:,:,oddExcl))+...
        1.005*squeeze(img1(:,:,:,evExcl));
    
    dw_mean(:,:,:)=mean(subtr,4);
% end
disp('++')
disp(['++ Size: ' num2str(size(dw_mean)) ' ++'])
disp('++')
%%
figure
plot((real((fft(squeeze(conj(dw_mean(19,22,:))))))));
saveas(gcf,['Exlc_' outSuffix '_example_19_22.png'],'png')
close gcf
%%
clear dw_mean_hlsvd
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
GF=10000;
LB=0;
disp('++ Begin HLSVD. This will take a while...')
for i=1:48
    for k=1:48
[dw_mean_hlsvd(i,k,:,c),~,~,~,~,~] = HLSVDPROrec((squeeze(dw_mean(i,k,:,c))),M,fsampl,t/1000,MM,freqrange);
    end
end
%end
%%
figure
plot((real((fft((squeeze(conj(dw_mean_hlsvd(22,19,:)))))))));
saveas(gcf,['Excl_' outSuffix 'post_hlsvd_19_22.png'],'png')
close gcf
%% still needs image water to be updated.

%%
save([subj outSuffix '_Excl.mat'],'dw_mean_hlsvd','img_water','meta_mask','mvmt','img1','-mat');
%end
    genspa(subj,[subj outSuffix '_Excl.mat'],[outSuffix '_Excl'])
