function  plotMRSIoverlay(metab_map2,T1_img,slice,metabnum,clip,fig)
%%

%% colormap
temp = metab_map2(:,:,slice,metabnum,1);
temp = temp(:);
if isempty(clip) == 0
    temp(temp < min(clip)) = min(clip);
    temp(temp>max(clip)) = max(clip);
end
temp = floor(temp*1000);

clrmap = jet(max(max(temp)));

im = metab_map2(:,:,slice,metabnum,1);
%%
imalpha = repmat(0.6,size(im,1),size(im,2));
imalpha(im == 0) = 0;
if isempty(clip) == 0
    im(im < min(clip)) = min(clip);
    im(im>max(clip)) = max(clip);
end
disp(['Image max = ' max(max(im))])
im = floor(im*1000);
%%


%% create look up table for indexing


%% resizing and colorizing
imout = ind2rgb(im,clrmap);

%%

%% Make new figure if requested
if fig == 1
    figure('position',[0 0 1000 2000])
    axis off
end

%% plot T1
imagesc(T1_img(:,:,slice))
axis square
colormap('grAy')
hold on

%% plot metab
a = imagesc(imout);
set(a,'alphadata',double(imalpha))
%% clean up
view([90 -90])

%% display some img parameters
disp(['++ Image min = ',num2str(min(min(im(im>0)))/1000)])
disp(['++ Image max = ',num2str(max(max(im))/1000)]);
%% clean up for output
% mapParams = struct;
% mapParams.imScale = imScale;
% mapParams.img = imout;
