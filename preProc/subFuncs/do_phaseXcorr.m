function [img1] = do_phaseXcorr(abc_dw,LB,GF)
%%
%%
% wrapper for Eddy current correction and phase cross correlation.
%%
dw=1/1250;
tt=0:dw:dw*(1024-1);
LB1=LB;
numImgs = size(abc_dw,4);
img1 = zeros([size(abc_dw,1) size(abc_dw,2) size(abc_dw,3)+512 numImgs]);
aa = zeros(size(abc_dw,3)+512,numImgs); %create variable for xcorr and make room for zero padding
for kk=1:size(abc_dw,1)
    for jj=1:size(abc_dw,2)
        for zz=1:numImgs
            if (mod(zz,2))
                cc1=[conj(squeeze(abc_dw(kk,jj,:,zz))); zeros(512,1)]; %zero padding
                aa(:,zz)=cc1.*(exp(-tt*pi*LB-tt.^2/(GF^2)))';
            else
                cc1=[conj(squeeze(abc_dw(kk,jj,:,zz))); zeros(512,1)]; %zero padding

                aa(:,zz)=cc1.*(exp(-tt*pi*LB1-tt.^2/(GF^2)))';
            end
         end
         bb = freq_xcorr_uzay(aa ,1250,123.2);%works better

                  bb=phase_max_uzay(bb);
          img1(kk,jj,:,:) =  bb;      
    end

end