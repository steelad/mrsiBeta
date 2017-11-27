function [roiValsOut,roiValsMean] = applyROI_MRSI(Map,LUT,ROI,filename);
%% Apply ROI to new map %%
% inputs:
% - Map:            Map to take values from (can be map of any type, incl. T1 or
% segmented image
% - LUT:            Look up table
% - ROI:            ROI that you want to apply (should be a 3dmatrix
% - filename:       filename for writing out the ROI values.
% 
% outputs:
% - roiValsOut:     Values from each voxel
% - roiValsMean:    Mean values from each voxel
%
%

%%
temp = unique(LUT(ROI == 1));
vxlInds = temp(temp ~= 0);
roiValsOut = zeros(length(vxlInds),size(Map,4)+1);
disp(['++ Outsize = ' num2str(size(roiValsOut))]);
clear temp
                                                
                                                roiValsOut(:,1) = vxlInds;
for i = 1:size(Map,4)
    for j = 1:length(vxlInds)
        temp = Map(:,:,:,i);
        vals = temp(LUT == vxlInds(j));
        roiValsOut(j,i+1) = unique(vals);
    end
end
                                                
roiValsMean = mean(roiValsOut(:,2:end),1);
if ~isempty(filename)
    dlmwrite([filename '_voxelwise.csv'],metabVals,'delimiter',',','precision',3)
    dlmwrite([filename '_avg.csv'],metabValMean,'delimiter',',','precision',3)
    disp(['++ Wrote ROI data to the file ' filename])
    disp(['++ Wrote mean of ROI data to single-line file ' filename])
end
        
