function plotMRSIoverlay_Leg(metab_map2,T1_img_slice,metabnum,clip)
%% make square image because it is simpler for now
% T1img = subj1(:,:,140);
% T1img = rot90(T1img);
%%
insize = size(T1_img_slice);
finalSqSize = max(insize);
BigDim = insize == finalSqSize;
%% Add to make the sides even numbers
if mod(finalSqSize,2) == 1
    disp('++ Largest dimension not even. adding additional zeros-vector before padding.')
    finalSqSize = finalSqSize + 1;
    padding = zeros(1,min(insize));
    if find(BigDim == 0) == 1
        disp('1')
        T1_img_slice = horzcat(T1_img_slice,rot90(padding));
    else
        disp('2')
        T1_img_slice = vertcat(T1_img_slice,padding);
    end
insize = size(T1_img_slice); 
disp(['New insize = ',num2str(insize)]);
end
%% pad the image
sizeDiffT1 = finalSqSize - min(insize);
%%
padding = zeros(sizeDiffT1/2,finalSqSize);
%%
if find(BigDim == 0) == 1
    T1_img_slice = vertcat(padding,T1_img_slice,padding);
else
    padding = rot90(padding);
    T1_img_slice = horzcat(padding,T1_img_slice,padding);
end
%% clean up
clear padding sizeDiffT1 insize
%%
%% Set up metabolite map
%%
temp = metab_map2(:,:,metabnum,1);
temp = temp(:);
if isempty(clip) == 0
    temp(temp < min(clip)) = min(clip);
    temp(temp>max(clip)) = max(clip);
end
temp(temp == 999) = 0;
temp = floor(temp*1000);
%% to make hot colormap %%
clrmap = jet(max(max(temp))*2);
clrmap = clrmap(length(clrmap)/2:length(clrmap),:);

im = metab_map2(:,:,metabnum,1);
%% circular shift for tumor patient.
im = circshift(im,3,1);
im = circshift(im,-1,2);
%%
im(im == 999) = 0;
if isempty(clip) == 0
    im(im < min(clip)) = min(clip);
    im(im>max(clip)) = max(clip);
end
disp(['Image max = ' num2str(max(max(im)))])
im = floor(im*1000);


%% resizing and colorizing
imout = ind2rgb(im,clrmap);
imScale = finalSqSize/size(im,1);
imout = imresize(imout,imScale,'nearest');
imalpha = repmat(0.6,size(im,1),size(im,2));
%%
imalpha(im == 0) = 0;
imalpha = imresize(imalpha,imScale,'nearest');

imagesc(T1_img_slice)
axis square
colormap('grAy')
hold on
%%
a = image(imout);
%%
set(a,'alphadata',double(imalpha))
