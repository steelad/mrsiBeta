function [LUTo, T1, Map1, Map2] = setUpMRSIWorkspace_LUT...
    (mapFile,T1file,LUTfile,trans);
%% Sets up workspace for MRSI if you have reconstructed your slab from a .rda file.
% inputs:
%
%- mapfile = metab_map workspace output from visualize LCModel scripts.
% assumes metab_map1 and metab_map2 are the names of the maps that should
% be created.
% - T1 file - nifti file for T1. Will have the orietation information
% stripped
% - LUTfile - reconstructed slab. Will have the orientation information
% stripped
%
% outputs:
%
% - LUTo - Lookup table for MRSI functions
% - T1 - T1 file for MRSI functions
% - Map1 & Map2 are the map1 and map2 files, respectively.
%%
load(mapFile)
ready = input('++ Enter 1 to remove orientation from the LUT and T1:  ');
if ready == 1
                cmd1 = ['fslorient -deleteorient ' T1file];
                disp(cmd1)
                system(cmd1); 
                disp('++ Orientation information removed from T1.')
                cmd2 = ['fslorient -deleteorient ' LUTfile];
                disp(cmd2)
                system(cmd2);
                disp('++')
                disp('++ Orientation information removed from LUT.')
                
            
end
    global AF_PATH
try
    isstring(AF_PATH);
catch
    AF_PATH = '';
    setenv('AF_PATH',AF_PATH)
end

if trans == 1
    ready = input('++ Enter 1 to confirm transpose metabolite maps:  ');
    if ready == 1
        trans == 1;
    else
        return
    end
end


setenv('LUT',LUTfile)
setenv('T1',T1file)
!echo $LUT
!${AF_PATH}3dcalc -a  $LUT  -expr 'a' -prefix LUT    
!${AF_PATH}3dcalc -a  $T1  -expr 'a' -prefix T1    
LUT = struct;
[~,LUT.img,LUT.info,~] = BrikLoad('LUT+orig');
[~,T1,~,~] = BrikLoad('T1+orig');
LUTo = LUT.img; 
    %%
    %%
    %%
    [~] = writeFilledLUT(LUT,metab_map1,...
        'filledMap1',trans);
    [~,Map1,~,~] = BrikLoad('filledMap1+orig');
    %Map1 = Map1.img;
    %%
    %%
    [~] = writeFilledLUT(LUT,metab_map2,...
        'filledMap2',trans);
    [~,Map2,~,~] = BrikLoad('filledMap2+orig');
    %Map2 = Map2.img;
    %%
end