function outputImg = writeFilledLUT...
    (LUT,metab_map2,outputName,trans);
disp('++ Loading image')
% setenv('AF_PATH','/Applications/afni/')
% setenv('LUT',LUTnii)
% setenv('outputName',outputName)
% !echo $LUT
% !$AF_PATH/3dcalc -a  $LUT  -expr 'a' -prefix LUT
%LUT = load_nii(LUTnii);
%%
%[~,LUT,Info,] = BrikLoad('LUT+orig');
%%
outputImg = LUT.img;
out = LUT.info;
disp('++ Beginning image fill')
%%
temp = fillImg(LUT.img,metab_map2,trans);
out.RootName = [outputName '+orig'];
%%
dd = zeros (1,5);
dd(1:length(size(temp))) = size(temp);
out.DATASET_DIMENSIONS = dd; 

%%
opt = struct;
opt.Prefix = outputName;
setenv('outputName',outputName)
%%
WriteBrik(temp,out,opt);
%%
!echo ++ Writing $outputName+orig as a Nifti file.
!${AF_PATH}3dcalc -a $outputName+orig -expr 'a' -prefix $outputName.nii

