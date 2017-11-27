function [out] = LCModel_run_array_csi_invivo_slaser(subj,hlsvdFile,fileSTEM,TE,recon)

%dirs=uipickfiles('out','struct')
dirs = dir(fileSTEM);
disp(['++ Number of files = ' num2str(length(dirs))]);
curdir=pwd;
for i=1:size(dirs,1)
[~,NAME,~] = fileparts(dirs(i).name);
cd(curdir);
    if TE < 50
        runlcmodel_arrayed_csi_spa_invivo_slaser_TE30(NAME,TE);
    else 
        warning('++ Using long TE')
        runlcmodel_arrayed_csi_spa_invivo_slaser_TE110(NAME,TE);
    end
cd(curdir)
end

out = 1;
%%
if recon == 1
  
    visualize_LCModel(subj,...
        hlsvdFile,curdir)
end
disp(['++ Done ++'])
