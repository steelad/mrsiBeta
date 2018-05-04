function genspa(subj,filename,outSuffix)
%subj = 'HV_102' ;
%%
load('semi_laser_30_33.spa','-mat');
%rmpath('/home/fs0/asteel/scratch/uzayScripts/scripts_templates/scripts/')
%cd (strcat('~/scratch/mrsiData/',subj))
%%

    load(filename,'-mat')
system(['mkdir ' outSuffix])
try 
    size(meta_mask);
catch
    load([subj '_finalWksp512.mat'],'meta_mask')
end
%%
numImgs = size(dw_mean_hlsvd,4);
%%

%% Generate spa file for each voxel in each image

            img_water = rot90(img_water);
 for MM = 1:numImgs
     disp(['++ Starting ' num2str(MM)])
 for i=1:48
     for k=1:48
         if(meta_mask(i,k)==1)
            aa=(conj(squeeze(dw_mean_hlsvd(i,k,:,MM))));
            aa = rot90(aa);
            data.metab=(conj((aa(1:512)')));
 
         
         data.ws=conj(squeeze(img_water(i,k,1:512,MM)));
         aa=ifft(fftshift(fft(data.ws)));
         data.ws=aa;
         B=[ subj '_' int2str(i) '_' int2str(k) '_' num2str(MM) '_' outSuffix '.spa'];
         save([outSuffix '/' B], 'data', '-mat')
       
         end
         
    end
 end
 
 end
 
disp('++ Generated all spa files. Switch to LCModel computer to fit.')
 %% save the final workspace
%%
