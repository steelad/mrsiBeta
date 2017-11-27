function resampleFull(inputFileStem,T1filename,outPrefix)

% inputfilename = 'test.nii'; T1filename = 'images_002_t1mprax1mmisowithNose32ch1001.nii.gz';
% outname = 'test_re.nii';
% fullCmd_do = ['3dWarp -deoblique -prefix temp ' inputfilename ];
% system(fullCmd_do)
% fullCmd = ['3dresample -overwrite -prefix ' outname ' -master ' T1filename ...
%     ' -inset ' inputfilename ' -rmode ' resamp]
setenv('t1',T1filename)
disp(['++ BETing'])
!bet  $t1 T1_ss.nii.gz -m -f 0.4
%%
disp('++ Doing Map 2.')
Map2 = resample_spectroscopy_slice( 'T1_ss.nii.gz','T1_ss_mask.nii.gz',...
    [inputFileStem '_MAP2.nii'], [outPrefix '_Map2_re']);
disp('++ Doing Map 1.')
Map1 = resample_spectroscopy_slice( 'T1_ss.nii.gz','T1_ss_mask.nii.gz',...
    [inputFileStem '_MAP1.nii'], [outPrefix '_Map1_re']);
disp('++ Doing LUT.')
LUT = resample_spectroscopy_slice( 'T1_ss.nii.gz','T1_ss_mask.nii.gz',...
    [inputFileStem '_LUT.nii'], [outPrefix '_LUT_re']);
disp('++ Doing GG.')
GG = resample_spectroscopy_slice( 'T1_ss.nii.gz','T1_ss_mask.nii.gz',...
    [inputFileStem '_GG.nii'], [outPrefix '_GG_re']);

%%
disp(['++ Done ++'])
%%
% system(fullCmd)
