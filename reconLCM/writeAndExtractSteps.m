%% steps to get images into matlab for value extraction
T1filename = 'images_002_t1mprax1mmisowithNose32ch1001.nii.gz';

[T1,Map,LUT] = writeResampledNiftis(metab_map2,T1filename,'102','test2');
plotMRSIoverlay_beta(Map,T1,140,11,[0 10],1);
getVoxelMetabMRSI_beta(Map,T1,LUT,140,[11 8 12 13])
%%
writeFilledLUT('LUT.nii',metab_map2,'testout.nii')