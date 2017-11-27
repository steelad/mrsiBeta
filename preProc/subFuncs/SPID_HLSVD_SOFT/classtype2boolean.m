function beta = classtype2boolean(position,classtype,opt)
%function that creates a matrix beat of booleans 
% beta(i,j)=1 if voxel i is similar to voxel j, =0 otherwise
% voxel i is similar to voxel j if the two voxels are close to each other
% (see iscloseto for the definition of closeness) AND both voxels have same
% tissue type (classtype)

bc = issameclass(classtype);
bs = iscloseto(position,opt);
beta = bs.*bc;