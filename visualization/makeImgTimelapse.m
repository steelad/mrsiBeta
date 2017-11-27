function makeImgTimelapse(img3d,img2,mov)
figure
colormap gray
axis square
% stdimg = std(img3d(:)); meanimg = mean(img3d(:));
if  length(size(img3d)) == 4
        disp('++ Maybe be a complex image...')
        if isreal(img3d) == 0
            img3d = squeeze(sum(real(abs(img3d)),3));
            if isempty(img2) == 0
                img2 = squeeze(sum(real(abs(img2)),3));
            end
        end
end
colorRange = [0 max(img3d(:))];

for i = 1:size(img3d,3)
    if isempty(img2) == 0
        subplot(1,2,1)
    end
  
    imagesc(img3d(:,:,i))
    caxis(colorRange)
    title(['Frame ' num2str(i)])
    axis square
    if isempty(img2) == 0
            subplot(1,2,2)
            imagesc(img2(:,:,i))
            caxis(colorRange)
            title(['Frame ' num2str(i) ' img2'])
            axis square
    end
    drawnow
    if isempty(mov) == 0
        mv(i) = getframe(gcf);
    else
        pause(0.1)
    end
end

    if isempty(mov) == 0
        writerObj = VideoWriter(mov,'Motion JPEG AVI');
        writerObj.FrameRate = 15;
        open(writerObj)
        writeVideo(writerObj,mv)
        close(writerObj)
    end
