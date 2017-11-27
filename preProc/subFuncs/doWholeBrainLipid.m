function [dw_mean_hlsvd,preLipidHLSVD,meta_mask,meta_mask2,lipid_mask] =...
    doWholeBrainLipid(subj,dw_mean_hlsvd,img_water,outName,TE)
%%

    temp = sum(abs(img_water),3); 
    disp('Carrying on through Lipid extraction because this is a whole brain image.')
        [meta_mask,meta_mask2,lipid_mask] = drawLipidMask(temp,dw_mean_hlsvd);
    preLipidHLSVD = dw_mean_hlsvd;
    [dw_mean_hlsvd] = doLipidSuppression(preLipidHLSVD,lipid_mask,meta_mask,1e-3);
    %%
    outSuffix = outName;
    save(outName,'preLipidHLSVD','dw_mean_hlsvd','img_water','meta_mask','-mat');
    for i = 1:size(dw_mean_hlsvd,1);
        for j = 1:size(dw_mean_hlsvd,2); 
            for ii = 1:size(dw_mean_hlsvd,4);
                dw_mean_hlsvd(i,j,:,ii) = conj(squeeze(dw_mean_hlsvd(i,j,:,ii)));
            end
        end
    end
    
    if TE == 110
        genspa_TE110(subj,outName,outSuffix)
    else
        genspa(subj,outName)
    end
    
    
