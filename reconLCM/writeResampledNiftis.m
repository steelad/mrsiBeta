function [T1,Map2,Map1,GluGabaMap,LUT] = ...
    writeResampledNiftis(subj,metab_map,T1_dir,outputName,sel)
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
%setenv('T1',T1_ss)
%% make spm T1
if isempty(T1_dir)
    P2 = spm_select(inf,'any')
    hdr2 = spm_dicom_headers(P2);
    disp(['++ Writing T1 imaging'])
    k2 = spm_dicom_convert(hdr2,'all','flat','nii')
else
end
%%
disp('++')
disp('++ Cleaning up, making T1 folder')

setenv('tempT1',k2.files{1})
    !mkdir T1
    !mv $tempT1 T1/T1.nii
disp('++')
disp('++ Skull Stripping')
    !3dSkullStrip -prefix T1/T1_ss.nii -input T1/T1.nii
    
%%


if isempty(sel)
    P1 = spm_select(1,'any'); %select spectroscopy dicom
else
    P1 = sel;
    disp(P1);
end

hdr1 = spm_dicom_headers(P1);
[~] = MRSI_write_wrapper(hdr1,metab_map,'__tempMAP2.nii');
[~] = MRSI_write_wrapper(hdr1,metab_map1,'__tempMAP1.nii');
[~] = MRSI_write_wrapper(hdr1,gluGabaCorrMap,'__tempGG.nii');
[~] = MRSI_write_wrapper(hdr1,LUT,'__tempLUT.nii');
[~] = MRSI_write_wrapper(hdr1,B,'C.nii');

%%
!3dcalc -a T1/T1_ss.nii -expr 'a' -prefix T1/T1_ss
!3dcalc -a T1/T1_ss+orig -expr 'a' -prefix ./T1_ss
T1file = 'T1_ss+orig';

resampleFull('__tempMAP.nii',T1file,'NN',[subj '_' outputName '_Map2']);
resampleFull('__tempMAP1.nii',T1file,'NN',[subj '_' outputName '_Map1']);
resampleFull('__tempGG.nii',T1file,'NN',[subj '_' outputName '_GluGabaMap']);
resampleFull('__tempLUT.nii',T1file,'NN',[subj '_LUT']);
%%
LUT = struct;
[~,T1,~,~] = BrikLoad(T1file);
[~,Map2,~,~] = BrikLoad([subj '_' outputName '_Map2+orig']);
[~,Map1,~,~] = BrikLoad([subj '_' outputName '_Map1+orig']);
[~,GluGabaMap,~,~] = BrikLoad([subj '_' outputName '_GluGabaMap+orig']);
[~,LUT.img,LUT.info,~] = BrikLoad([subj '_LUT+orig']);
%%
%%
save([subj '_' outputName '_workspace.mat'],'-v7.3')
delete __tempMAP.nii __tempLUT.nii __tempMAP.mat __tempLUT.mat
delete __tempMAP1.nii __tempGG.nii __tempMAP1.mat __tempGG.mat
