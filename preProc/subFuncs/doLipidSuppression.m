function [postLipidImg] = doLipidSuppression(dw_mean_hlsvd,lipid_mask,meta_mask,lipidPenalty)
disp('++ Setting up lipid script')
for i=1:48
    for j=1:48
        for ii = 1:size(dw_mean_hlsvd,4)
        temp(i,j,:,ii)=fftshift(fft(conj(squeeze(dw_mean_hlsvd(i,j,:,ii))))); % FFT to get into frequency domain, 
                                                                        % necessary for lipid 
        
        end
    end
end
img_metab_fft=temp;
N = size(img_metab_fft);
csi_rad = 48;                     % size of high res image, has lipid ringing
mask_csi = repmat(circular_mask([48,48],csi_rad),[1,1,N(3)]);
N
%% ste up parameters for the lipid suppression
%imagesc(mask_csi(:,:,1))
FT_csi = FT_v2(mask_csi);
img_zf = FT_csi'*(FT_csi*img_metab_fft);
x = img_zf;
%%
plot_projection([(img_zf), (img_zf)],2)
param = init;
param.xfmWeight = lipidPenalty;           % Lipid basis penalty
param.FT = FT_csi;
param.Itnlim = 20;
param.data =  param.FT*img_metab_fft;
param.Itnlim = 10;
param.Lipid = get_LipidBasis(img_zf,lipid_mask);    % Lipid basis functions obtained from low-res image
param.Bmask = meta_mask;                            % Brain mask
plot_projection([img_zf,img_metab_fft],1)
%%
%%
disp('++ Doing lipid suppression. This will take a while.')
disp('... Seriously. Go grab a coffee or something.')
tic
iters = 10
for t = 1:iters
    disp(['Starting iteration ', num2str(t), ' of ', num2str(iters)])
    x = lipid_suppression(x,param);
    figure(1), imagesc( [ meta_mask.*sum(abs(x),3) , meta_mask.*sum(abs(img_zf),3) ] ),...
    axis image
    colorbar
    drawnow
end
toc
%%
x_meta = x .* repmat(meta_mask,[1,1,N(3)]);

%overplot_spectra( img_zf .* repmat(meta_mask,[1,1,N(3)]), x_meta, 18, 17, 3, 3,512, 1024)
%%
clear temp

for i=1:48
    for j=1:48
        for ii = 1:size(dw_mean_hlsvd,4);
            img_metab_dw_90_fft_lipid(i,j,:,ii)=squeeze(x_meta(i,j,:,ii));
            
            temp(i,j,:,ii)=(ifft(fftshift(img_metab_dw_90_fft_lipid(i,j,:,ii))));
        end
    end
    
    
end
%%
postLipidImg = temp;
