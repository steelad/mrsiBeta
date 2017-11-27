function [metabVals,metabValMean,roiMask] = ...
    getVoxelMetabMRSI(metab_map_2,T1_img,LUT,...
                      slice,metabnum,filename);
%% Select voxels and draw regions on interest, extract metabolite data from input metabolite map
% inputs:
%  - metab_map_2 = metabolite map to select from
%  - T1_img = T1 image as underlay
%  - LUT = Look up table
%  - slice = slice number to draw (assumed axial orietntation for now)
%  - metab num = list of metabolite numbers (1..n).  First in list is used as overlay for drawing roi/
%  - filename = write out files with the name here
% outputs:
%  - metabVals = values from each voxel for all metabs (num vxl x num
%  metabs). Also written to file if desired.
%  - metabValMean = mean values from roi for all metabs (1 x num metabs).
%  Also written to file if desired
%  - roiMask = ROI mask to apply to other images. Also written to file if
%  desired

%%
%%
%% detect clip
metabSlice = metab_map_2(:,:,slice,metabnum(1));
T1_img_slice = T1_img(:,:,slice);
temp = metabSlice(metabSlice > 0);
mT = mean(temp);
stdT = std(temp);
uClip = mT+2*stdT;
%% Get indices of interest with user input
plotMRSIoverlay(metabSlice,T1_img_slice,1,1,[0 uClip],1);
axis square
drawnow
[Vx,Vy] = ginput; 

%% translate to MRSI image
vxlInds = zeros(length(Vx));
Vx = round(Vx); Vy = round(Vy);
for i = 1:length(vxlInds)
    vxlInds(i) = LUT(Vy(i),Vx(i),slice);
end
vxlInds = unique(vxlInds);
vxlInds = vxlInds(vxlInds ~= 0);
%% get metabolite values
mask = zeros(size(metabSlice));
metabVals = zeros(length(vxlInds),size(metab_map_2,4)+1);
    roiMask = zeros(size(LUT));
metabVals(:,1) = vxlInds;
for i = 1:length(vxlInds)
    mask(LUT(:,:,slice) == vxlInds(i)) = 1;
    for j = 1:size(metab_map_2,4)
        temp = metab_map_2(:,:,slice,j);
        temp1 = temp(LUT(:,:,slice) == vxlInds(i));
        temp1 = unique(temp1);
        metabVals(i,j+1) = temp1;
    end
    roiMask (LUT == vxlInds(i)) = 1;
         
end
% for i = 1:size(metab_map_2,4)
%     temp = 
%     temp1 = temp(mask == 1);
%     temp1 = unique(temp1);
%     temp1(temp1~=0);
%     clear temp1
% end
% metabVals = unique(metabVals);
metabValMean = mean(metabVals(:,2:end),1);
%% 
%% plot result
close gcf
%% Plot Metab
clrs = lines(length(metabnum));
figure('position',[300 500 150*length(metabnum) 600])
hold on
for i = 1:length(metabnum)
plot([i-.25 i+.25],[metabValMean(metabnum(i)) metabValMean(metabnum(i))],'-','linewidth',2,'color',clrs(i,:));
scatter(repmat(i,1,size(metabVals,1)),metabVals(:,metabnum(i)+1),50,'linewidth',2,'markeredgecolor',clrs(i,:));
end
%xlim([0.5 1.5])
set(gca,'linewidth',2,'xtick',1:length(metabnum),'xticklabel',metabnum,'fontsize',16)
%%

figure('position',[300+150*length(metabnum) 500 600 600])
plotROIslice(T1_img_slice,mask,0)
%%
if isempty(filename) ~= 1
    opt.Prefix = [filename '_roi'];
    [~,~,info,~] = BrikLoad('T1_ss+orig');
    WriteBrik(roiMask,info,opt)
    setenv('outputName',filename)
    !${AF_PATH}3dcalc -a  ${outputName}_roi+orig  -expr 'a' -prefix ${outputName}_roi.nii
    dlmwrite([filename '_voxelwise.csv'],metabVals,'delimiter',',','precision',3)
    dlmwrite([filename '_avg.csv'],metabValMean,'delimiter',',','precision',3)
    disp(['++ Wrote ROI data to the file ' filename])
    disp(['++ Wrote mean of ROI data to single-line file ' filename])
    
end


