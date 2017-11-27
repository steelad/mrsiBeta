function plotROIslice(T1_img_slice,mask,fig)
%% plot T1 
if fig == 1
    figure('position',[500 500 600 600])
end
imagesc(T1_img_slice)
axis square off
colormap('grAy')
hold on

%% plot metab
a = image(cat(3,mask,zeros(size(mask)),zeros(size(mask))));
set(a,'alphadata',double(mask))
%%
view([90 -90])
