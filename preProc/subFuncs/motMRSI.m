function [img_dw_90_water,img_dw_90,mvmt] = motMRSI(imgWater,img_dw_90_water,img_dw_90,mot)
    if mot == 1
        interpMet = 'nearest';
    elseif mot == 2
        interpMet = 'linear';
    end
    numImgs = size(imgWater,4);
    disp('++ Start motion correction with rigid body registration')
    %% register images
        transMat = cell(1,numImgs);
        [optim,met] = imregconfig('monomodal');
        prop = imref2d(size(imgWater(:,:,1)));
        water_img1 = zeros(size(imgWater));
        mvmt = zeros(3,3,numImgs);
    for i = 1:numImgs
        displacementOutmatrix = imregtform(imgWater(:,:,i),imgWater(:,:,1),...
            'rigid',optim,met);
        transMat{i} = displacementOutmatrix;
        [water_img1(:,:,i),~] = imwarp(imgWater(:,:,i),...
            prop,displacementOutmatrix,'outputview',prop,'interp',interpMet);
        mvmt(:,:,i) = transMat{i}.T;
    end
%%
mvmt = reshape(mvmt,[],10);
figure
subplot(2,1,1); title(['rotation (q)'])
plot(sin(mvmt(2,:)),'b')
subplot(2,1,2); title(['translation (total, vxl)'])
plot(mvmt(3,:).^2.*mvmt(6,:).^2,'r')
saveas(gcf,['motion_' outName '.png'],'png')
close gcf
imgWater = img_dw_90_water;
for i = 1:numImgs
    disp(['++ Applying registration: ' num2str(i)])
    for j = 1:numTPs
        img_dw_90(:,:,j,i) = imwarp(img_dw_90(:,:,j,i),...
        prop,transMat{i},'outputview',prop,'interp',interpMet);
        img_dw_90_water(:,:,j,i) = imwarp(img_dw_90_water(:,:,j,i),...
        prop,transMat{i},'outputview',prop,'interp',interpMet);
    end
end
makeImgTimelapse(imgWater,img_dw_90_water,'mvmt.mp4')
close all