function writeNiftis(subj,metab_map,outputName,sel)
%%
% This program takes the metabolite maps file and writes out metab_map2 and look up table.
% Should also include a wrapper to write read and write the rilled look up
% table with metab_map1, gluGabaCorrMap, etc.
%
%%
%%
load(metab_map)
metab_map = metab_map2;
LUT = 1:size(metab_map,1)*size(metab_map,2);
LUT = reshape(LUT,size(metab_map,1),size(metab_map,2));

if isempty(sel)
    P1 = spm_select(1,'any'); %select spectroscopy dicom
else
    P1 = sel;
    disp(P1);
end

hdr1 = spm_dicom_headers(P1);
[~] = MRSI_write_wrapper(hdr1,metab_map2,[subj '_' outputName '_MAP2.nii']);
[~] = MRSI_write_wrapper(hdr1,metab_map1,[subj '_' outputName '_MAP1.nii']);
[~] = MRSI_write_wrapper(hdr1,gluGabaCorrMap,[subj '_' outputName '_GG.nii']);
[~] = MRSI_write_wrapper(hdr1,LUT,[subj '_' outputName '_LUT.nii']);
%%

delete *_MAP.mat *_LUT.mat
delete *_MAP1.mat *_GG.mat
